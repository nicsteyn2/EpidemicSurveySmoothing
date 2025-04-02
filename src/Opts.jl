
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
    end

    return(opts)

end