
using CSV, DataFrames

include("../src/Simulations.jl")
include("../src/Opts.jl")
include("../src/Models.jl")
include("../src/FitModel.jl")
include("../src/LoadData.jl")

start_round = 1
end_round = 19
n = 5 # Target successful runs

fname = "paper/outputs/s8-timing/steynRt_R$(start_round)to$(end_round)_2000.csv"

# Load data
(Y, tMax) = loadReactData(start_round, end_round)

# Load PCR curve
pcr_raw = CSV.read("data/PCR-curves/hellewell.csv", DataFrame)
pcr = vec(mean(Matrix(pcr_raw[:,2:end]), dims=1))
pcr_full = zeros(tMax)
pcr_full[1:min(length(pcr), tMax)] = pcr[1:min(length(pcr), tMax)]

# Specify generation time distribution
gentime = pdf.(Gamma(2.29, 1/0.36), 0:tMax) 
gentime = gentime/sum(gentime)

# Specify the initial distribution for Pt
phat = max(Y.nPos[1], 1)/Y.nSamples[1]
sd = sqrt((phat*(1-phat))/Y.nSamples[1])
pP0 = Truncated(Normal(phat, sd), 0, 1.0)

# Set model options
opts = fetchDefaultOpts("RtModel")
opts["N"] = 2000
opts["L"] = 30
opts["L-marginalposterior"] = min(50, tMax)
opts["maxRhat"] = 1.05
opts["minESS"] = 100
opts["pP0"] = pP0
opts["pR0"] = Uniform(0.4, 1.2)
opts["genTimeDist"] = gentime
opts["infToPositivity"] = pcr_full





# Function to fit the model
function fit_model_once(Y, opts)

    # Fit the model
    st = time()
    (df_states, df_params, θ, diag) = fitModel(RtModelForPaper, Y, deepcopy(opts); skipResamplingPMMH=false, verbose=false, posteriorPredictive="betabinom", checkLogLikStdDev=false)
    en = time()

    maxrhat = maximum(df_params.rhat)
    miness = minimum(df_params.ess)
    time_taken = en - st

    return(maxrhat, miness, time_taken)

end


function fit_model_ntimes(Y, opts, n, fname)

    # Pre-allocate output dataframe
    df_out = DataFrame()

    n_attempts = 1
    n_successful = 0

    while ((n_attempts <= 50) & (n_successful < n))

        # Fit the model
        (maxrhat, miness, time_taken) = fit_model_once(Y, opts)

        # Check if successful
        is_success = (maxrhat <= 1.05) & (miness >= 100)

        # Store outputs
        df_temp = DataFrame(start_round=Y.Round[1], end_round=Y.Round[end], iter=n_attempts, success=is_success, maxrhat=maxrhat, miness=miness, time_taken=time_taken, num_obs = sum(Y.nSamples .> 0), num_days = nrow(Y))
        df_out = vcat(df_out, df_temp) 
        CSV.write(fname, df_out)
        display(df_out)
                            
        n_attempts += 1
        n_successful += is_success

    end

end


fit_model_ntimes(Y, opts, n, fname)