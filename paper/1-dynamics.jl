
using CSV, DataFrames

include("../src/Simulations.jl")
include("../src/Opts.jl")
include("../src/Models.jl")
include("../src/FitModel.jl")

function generateConfigTable()
    
    # Set default options
    rho = [0, 0, 2e-4]
    xi = [0, 0.5, 0]
    nSamples = [10, 100, 1000, 10000, 100000]

    # Generate config table for varying nt (low sigma)
    df_config_low = DataFrame()
    df_config_low[!,:nSamples] = repeat(nSamples, inner=3)
    df_config_low[!,:rho] = repeat(rho, outer=length(nSamples))
    df_config_low[!,:xi] = repeat(xi, outer=length(nSamples))
    df_config_low[!,:sigma] .= 0.008
    
    # Generate config table for varying nt (high sigma)
    df_config_high = DataFrame()
    df_config_high[!,:nSamples] = repeat(nSamples, inner=3)
    df_config_high[!,:rho] = repeat(rho, outer=length(nSamples))
    df_config_high[!,:xi] = repeat(xi, outer=length(nSamples))
    df_config_high[!,:sigma] .= 0.016
    
    # Combine the tables and ignore duplicate rows
    df_config = unique(vcat(df_config_low, df_config_high))
    df_config[!,:confignum] = 1:nrow(df_config)
    
    # Save the configuration table
    CSV.write("paper/1-dynamics-configtable.csv", df_config)
    
end

function simulateData(nSamples, sigma, rho, xi)
    
    # Simulate 100 days of an epidemic
    (rt, Pt) = simulateBoundedRandomWalk(tMax=100, σ=sigma)
    
    # Simulate a survey
    if (rho > 0) & (xi > 0)
        error("Cannot simulate from the beta-binomial and weighted models simultaneously.")
    elseif rho > 0
        Y = simulateSurvey(nSamples, rt, Pt; method="betabinom", params=[rho])
    elseif xi > 0
        Y = simulateSurvey(nSamples, rt, Pt; method="weighted", params=[xi])
    else
        Y = simulateSurvey(nSamples, rt, Pt; method="binomial")
    end
    
    return(Y)
    
end


function fitModels(Y, sigma, rho, xi)
    
    # Define output dataframes
    df_states_out = DataFrame()
    df_params_out = DataFrame()
    
    # Specify the initial distribution for Pt
    phat = max(Y.nPos[1], 1)/Y.nSamples[1]
    sd = sqrt((phat*(1-phat))/Y.nSamples[1])
    pP0 = Truncated(Normal(phat, sd), 0.001, 0.3)
    
    # Fit binomial model at unknown parameters
    println("Starting binomial model (unknown parameters)...")
    opts = fetchDefaultOpts("binomial")
    opts["pP0"] = pP0
    opts["posteriorNumberOfParticles"] = 10000
    st = time()
    (df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=false, posteriorPredictive="binomial")
    en = time()
    df_states[!,:model] .= "binomial"
    df_states[!,:params] .= "estimated"
    df_params[!,:model] .= "binomial"
    df_params[!,:params] .= "estimated"
    df_params[!,:runtime] .= en - st
    df_states_out = vcat(df_states_out, df_states)
    df_params_out = vcat(df_params_out, df_params)
    
    # Fit binomial model at known parameters
    println("Starting binomial model (known parameters)...")
    opts["N"] = 100000
    df_states = fitModelFixedParams([sigma], BinomialModel, Y, opts; posteriorPredictive="binomial")
    df_states[!,:model] .= "binomial"
    df_states[!,:params] .= "fixed"
    df_states_out = vcat(df_states_out, df_states)
    
    # Fit beta-binomial model at unknown parameters
    println("Starting beta-binomial model (unknown parameters)...")
    opts = fetchDefaultOpts("betabinom")
    opts["pP0"] = pP0
    opts["posteriorNumberOfParticles"] = 10000
    st = time()
    (df_states, df_params, θ, diag) = fitModel(BetaBinomialModel, Y, opts; skipResamplingPMMH=true, verbose=false, posteriorPredictive="betabinom")
    en = time()
    df_states[!,:model] .= "betabinom"
    df_states[!,:params] .= "estimated"
    df_params[!,:model] .= "betabinom"
    df_params[!,:params] .= "estimated"
    df_params[!,:runtime] .= en - st
    df_states_out = vcat(df_states_out, df_states)
    df_params_out = vcat(df_params_out, df_params)
    
    # Fit beta-binomial model at known parameters
    println("Starting beta-binomial model (known parameters)...")
    opts["N"] = 100000
    if rho > 0
        df_states = fitModelFixedParams([sigma, rho], BetaBinomialModel, Y, opts; posteriorPredictive="betabinom")
    else
        df_states = fitModelFixedParams([sigma], BinomialModel, Y, opts; posteriorPredictive="binomial")
    end
    df_states[!,:model] .= "betabinom"
    df_states[!,:params] .= "fixed"
    df_states_out = vcat(df_states_out, df_states)
    
    # Fit weighted model at unknown parameters
    println("Starting weighted model (unknown parameters)...")
    opts = fetchDefaultOpts("weighted")
    opts["pP0"] = pP0
    opts["posteriorNumberOfParticles"] = 10000
    opts["boundPrevalence"] = true
    opts["minPrevalence"] = 0.001
    opts["maxPrevalence"] = 0.3
    st = time()
    (df_states, df_params, θ, diag) = fitModel(WeightedModel, Y, opts; skipResamplingPMMH=true, verbose=false, posteriorPredictive="weighted")
    en = time()
    df_states[!,:model] .= "weighted"
    df_states[!,:params] .= "estimated"
    df_params[!,:model] .= "weighted"
    df_params[!,:params] .= "estimated"
    df_params[!,:runtime] .= en - st
    df_states_out = vcat(df_states_out, df_states)
    df_params_out = vcat(df_params_out, df_params)
    
    # Note: we cannot fit the weighted model at known parameters
    
    return(df_states_out, df_params_out)
    
end

function fitMultipleModels(N, nSamples, sigma, rho, xi)
    
    df_states_out = DataFrame()
    df_params_out = DataFrame()
    df_data_out = DataFrame()
    
    ProgBar = Progress(N, dt=1, barlen=50, desc="Fitting $N models...")
    
    for ii = 1:N
        
        success = false
        counter = 0
        while !success && counter < 20
            try
                Y = simulateData(nSamples, sigma, rho, xi)
                (df_states, df_params) = fitModels(Y, sigma, rho, xi)
                df_states[!,:iter] .= ii
                df_params[!,:iter] .= ii
                Y[!,:iter] .= ii
                df_states_out = vcat(df_states_out, df_states)
                df_params_out = vcat(df_params_out, df_params)
                df_data_out = vcat(df_data_out, Y)
                next!(ProgBar)
                success = true
            catch e
                counter += 1
            end
        end
        
        if counter >= 20
            error("Failed to fit model $ii.")
        end
        
    end
    
    return(df_states_out, df_params_out, df_data_out)
    
end

function runConfigNtimes(N, confignum)
    
    # Load the configuration table
    config = CSV.read("paper/1-dynamics-configtable.csv", DataFrame)
    nSamples = config.nSamples[confignum]
    sigma = config.sigma[confignum]
    rho = config.rho[confignum]
    xi = config.xi[confignum]
    
    # Run the models
    (df_states, df_params, df_data) = fitMultipleModels(N, nSamples, sigma, rho, xi)
    
    # Append confignums
    df_states[!,:confignum] .= confignum
    df_params[!,:confignum] .= confignum
    df_data[!,:confignum] .= confignum
    
    # Save the output
    fname = "paper/outputs/1-dynamics/1-dynamics-config$confignum-"
    CSV.write(fname * "states.csv", df_states)
    CSV.write(fname * "params.csv", df_params)
    CSV.write(fname * "data.csv", df_data)
    
end
