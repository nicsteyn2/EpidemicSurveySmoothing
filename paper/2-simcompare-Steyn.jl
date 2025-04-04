
using ProgressMeter, CSV

include("../src/Simulations.jl")
include("../src/Opts.jl")
include("../src/Models.jl")
include("../src/FitModel.jl")


# Generate N simulations from the beta-binomial model
function saveSimulations(N; tMax=100, σ=0.01, ρ=2e-4, nSamples=5000)
    
    # Generate simulations
    Yall = DataFrame()
    ProgBar = Progress(N, dt=1, barlen=50, desc="Simulating...")
    for ii = 1:N
        (rt, Pt) = simulateBoundedRandomWalk(tMax=tMax, σ=σ)
        Y = simulateSurvey(nSamples, rt, Pt; method="betabinom", params=[ρ])
        Y[:,:iter] .= ii
        Yall = vcat(Yall, Y)
        next!(ProgBar)
    end
    
    # Tidy up formatting
    Y_growth = Yall[:, [:t, :rt, :iter]]
    rename!(Y_growth, Dict(:rt => :value))
    Y_growth[:,:variable] .= "rt"
    
    Y_prev = Yall[:, [:t, :Pt, :iter]]
    rename!(Y_prev, Dict(:Pt => :value))
    Y_prev[:,:variable] .= "Pt"
    
    Y_swabs = Yall[:, [:t, :nPos, :iter]]
    rename!(Y_swabs, Dict(:nPos => :value))
    Y_swabs[:,:variable] .= "nPos"

    Y_samples = Yall[:, [:t, :nSamples, :iter]]
    rename!(Y_samples, Dict(:nSamples => :value))
    Y_samples[:,:variable] .= "nSamples"
    
    Yout = vcat(Y_growth, Y_prev, Y_swabs, Y_samples)
    
    CSV.write("paper/outputs/2-simcompare/simulations_steyn.csv", Yout)
    return(Yout)
    
end



function loadSimulations()
    
    steyn = CSV.read("paper/outputs/2-simcompare/simulations_steyn.csv", DataFrame)
    abbott = CSV.read("paper/outputs/2-simcompare/simulations_abbott.csv", DataFrame)
    eales = CSV.read("paper/outputs/2-simcompare/simulations_eales.csv", DataFrame; missingstring="NA")
    replace!(eales.value, missing => NaN)
    
    steyn[:,:model] .= "Steyn"
    abbott[:,:model] .= "Abbott"
    eales[:,:model] .= "Eales"
    
    return(vcat(steyn, abbott, eales))
    
end



function fitAllSimulations(N)
    
    sims = loadSimulations()
    df_states_all = DataFrame()
    df_params_all = DataFrame()
    
    combinations = DataFrame(model=repeat(["Steyn", "Abbott", "Eales"], inner=N), iter=repeat(1:N, 3))
    
    ProgBar = Progress(nrow(combinations), dt=1, barlen=50, desc="Fitting models...")
    for ii in 1:nrow(combinations)
        
        # Set options
        iter = combinations.iter[ii]
        model = combinations.model[ii]
        println("Fitting model $model, iteration $iter")
        
        # Process into correct format
        Yin = sims[(sims.iter .== iter) .& (sims.model .== model), :]
        Y = DataFrame(t=1:100)
        Y.nSamples = Int.(Yin[Yin.variable .== "nSamples", :value])
        Y.nPos = Int.(Yin[Yin.variable .== "nPos", :value])
        Y.rt = Yin[Yin.variable .== "rt", :value]
        Y.Pt = Yin[Yin.variable .== "Pt", :value]
        
        # Specify the initial distribution for Pt
        phat = max(Y.nPos[1], 1)/Y.nSamples[1]
        sd = sqrt((phat*(1-phat))/Y.nSamples[1])
        pP0 = Truncated(Normal(phat, sd), 0, 1.0)
        
        # Load options dictionary
        opts = fetchDefaultOpts("betabinom")
        opts["simulatedData"] = false
        opts["showChunkProgress"] = false
        opts["pP0"] = pP0
        
        # Fit the model
        (df_states, df_params, θ, diag) = fitModel(BetaBinomialModel, Y, opts; skipResamplingPMMH=true, verbose=false, posteriorPredictive="betabinom", checkLogLikStdDev=false)
        
        # Save oututs
        df_states[:,:model] .= "Steyn"
        df_states[:,:simulation] .= model
        df_states[:,:iter] .= iter
        
        df_params[:,:model] .= "Steyn"
        df_params[:,:simulation] .= model
        df_params[:,:iter] .= iter
        
        df_states_all = vcat(df_states_all, df_states)
        df_params_all = vcat(df_params_all, df_params)

        next!(ProgBar)
        
    end
    
    CSV.write("paper/outputs/2-simcompare/simcompare_steyn_states.csv", df_states_all)
    CSV.write("paper/outputs/2-simcompare/simcompare_steyn_params.csv", df_params_all)
    
    return(df_states_all, df_params_all)
    
end

# # Generate simulated data
# sims = saveSimulations(100)

# Fit models to simulated data
(states_test, params_test) = fitAllSimulations(10)

