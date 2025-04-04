
using Distributions

function fetchDefaultOpts(model)
    
    opts = Dict(
    # SMC options
    "N" => 1000, # Number of particles to use
    "L" => 50, # Fixed-lag resampling length (when estimating the hidden states)
    "ω" => pdf.(Gamma(2.36, 2.74), 1:100), # Default generation time dist.
    "pr0" => Uniform(-0.3, 0.3), # Initial distribution for rt,
    "pP0" => missing, # Initial distribution for Pt (this depends on reported case data, specify at run-time),
    "stateNames" => ["rt", "Pt"],
    # PMMH options
    "nChains" => 3, # Number of PMMH chains to run
    "chunkSize" => 100, # Number of samples per chunk
    "maxChunks" => 50, # Maximum number of PMMH chunks to run
    "maxRhat" => 1.05,  # Stopping criterion: maximum Rhat value
    "minESS" => 100, # Stopping criterion: minimum effective sample size
    # Marginal posterior options
    "posteriorNumberOfParticles" => 2000,
    "posteriorParamSamples" => 100
    )
    
    if model == "binomial"
        # Set binomial-specific options
        # PMMH options
        opts["paramPriors"] = [Uniform(0, 0.2)]
        opts["initialParamSamplers"] = [Uniform(0, 0.2)]
        opts["paramLimits"] = [(0, 0.2)]
        opts["propStdDevInit"] = [0.01]
        opts["paramNames"] = ["sigma"]
    elseif model == "betabinom"
        # Set beta-binomial specific options
        opts["paramPriors"] = [Uniform(0, 0.2), Uniform(0, 0.01)]
        opts["initialParamSamplers"] = [Uniform(0, 0.2), Uniform(0, 1e-3)]
        opts["paramLimits"] = [(0, 0.2), (0, 0.01)]
        opts["propStdDevInit"] = [0.01, 5e-5]
        opts["paramNames"] = ["sigma", "rho"]
    elseif model == "weighted"
        # Set weighted model specific options
        opts["paramPriors"] = [Uniform(0, 0.2), Uniform(0.1, 10)]
        opts["initialParamSamplers"] = [Uniform(0.005, 0.015), Uniform(0.67, 1.5)]
        opts["paramLimits"] = [(0, 0.2), (0.1, 10)]
        opts["propStdDevInit"] = [0.005, 0.1]
        opts["paramNames"] = ["sigma", "c"]
        opts["boundPrevalence"] = false # Useful in low-data scenarios to stop the model getting "trapped" at low prevalence    end
    elseif model == "RtModel"
        opts["stateNames"] = ["Rt", "It", "Pt"]
        opts["paramPriors"] = [Uniform(0, 1), Uniform(0, 0.01)]
        opts["initPeriod"] = 10
        opts["initialParamSamplers"] = [Uniform(0.05, 0.2), Uniform(1e-4, 4e-4)]
        opts["paramLimits"] = [(0, 1), (0, 0.01)]
        opts["propStdDevInit"] = [0.04, 1e-4]
        opts["paramNames"] = ["sigma_R", "rho"]
    else
        error("Invalid model name.")
    end

    return(opts)
    
end