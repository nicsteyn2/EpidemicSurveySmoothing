
using CSV, DataFrames

include("../src/Simulations.jl")
include("../src/Opts.jl")
include("../src/Models.jl")
include("../src/FitModel.jl")
include("../src/LoadData.jl")

start_round = 17
end_round = 19

# Load data
(Y, tMax) = loadReactData(start_round, end_round)

# Specify the initial distribution for Pt
phat = max(Y.nPos[1], 1)/Y.nSamples[1]
sd = sqrt((phat*(1-phat))/Y.nSamples[1])
pP0 = Truncated(Normal(phat, sd), 0, 1.0)

# Set model options
opts = fetchDefaultOpts("binomial")
opts["showChunkProgress"] = false
opts["maxRhat"] = 1.05
opts["minESS"] = 100
opts["pP0"] = pP0

# Set appropriate N
if nrow(Y) >= 500
    # opts["propStdDevInit"] = [0.003, 5e-5]
    opts["posteriorNumberOfParticles"] = 2000
    opts["L-marginalposterior"] = 70
end

# # Check standard deviation of log lik manually
# (sd, ll) = estimateStdDevLogLik(100, WeightedModel, [0.02, 1.8], Y, opts)

# # Check model
# (X, W) = WeightedModel([0.015, 1], Y, opts)


# Fit the model
st = time()
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
en = time()

# Save outputs
CSV.write("paper/outputs/s9b-react-basic/steyn_r$(start_round)to$(end_round)_states.csv", df_states)
CSV.write("paper/outputs/s9b-react-basic/steyn_r$(start_round)to$(end_round)_params.csv", df_params)

# Calculate nPos coverage
l = df_states.lower[df_states.variable .== "nPos"]
u = df_states.upper[df_states.variable .== "nPos"]
inPredInterval = (Y.nPos .>= l) .& (Y.nPos .<= u)
coverage = mean(inPredInterval[Y.nSamples .> 0])

# Calculate width of credible intervals on Pt
l = df_states.lower[df_states.variable .== "Pt"]
u = df_states.upper[df_states.variable .== "Pt"]
width = u - l
meanwidth = mean(width[Y.nSamples .> 0])

# Save model summary
outfile = open("paper/outputs/s9b-react-basic/steyn_r$(start_round)to$(end_round)_summary.txt", "w")
write(outfile, "Summary of results...\n\n")
write(outfile, "Time taken = $(en-st) seconds\n\n")
write(outfile, "Max Rhat = $(round(maximum(df_params.rhat), digits=2))\n")
write(outfile, "Min ESS = $(round(minimum(df_params.ess), digits=2))\n\n")
write(outfile, "sigma = $(round(df_params.m[1], sigdigits=2)) ($(round(df_params.l[1], sigdigits=2)), $(round(df_params.u[1], sigdigits=2)))\n")
write(outfile, "Coverage of nPos = $(round(100*coverage, digits=1))%\n")
write(outfile, "Mean width of credible intervals on Pt (%) = $(round(100*meanwidth, digits=4))\n")
close(outfile)