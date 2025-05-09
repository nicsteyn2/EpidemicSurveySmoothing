
using MCMCChains, LinearAlgebra
include("Likelihood.jl")
include("Support.jl")


function PMMH(bootstrapFilter::Function, Y::DataFrame, opts::Dict; verbose=false)
    
    # Extract options
    nChains = opts["nChains"]
    chunkSize = opts["chunkSize"]
    maxChunks = opts["maxChunks"]
    nParams = length(opts["paramPriors"])
    propStdDev = opts["propStdDevInit"]

    maxSamples = maxChunks * chunkSize

    meanEst = zeros(nParams)'
    covEst = Diagonal(opts["propStdDevInit"].^2) 
    if "initCovEst" in keys(opts)
        covEst = opts["initCovEst"]
        opts["initCovWeight"] = 2000
    else
        opts["initCovWeight"] = 100
    end


    # Check that we are actually running multiple
    if nChains < 2
        error("These methods require a minimum of two PMMH chains to be run.")
    end
    
    # Pre-allocate output
    θall = zeros(maxSamples, nParams, nChains) # A nSamples x nParams x nChains array
    diag = DataFrame() # A diagnostics dataframe, note the different format to simplePMMH()
    
    (ii, ESS, Rhat, maxTuningInd, nTuningChunks, isTuning) = (1, 0, Inf, Inf, 0, true)
    while (ii <= opts["maxChunks"] && (Rhat > opts["maxRhat"] || ESS < opts["minESS"]))
        
        # Fetch inds to save to
        firstInd = (ii-1)*chunkSize + 1
        lastInd = ii*chunkSize
        inds = firstInd:lastInd
        
        # Run a single chunk of PMMH for each chain
        nAcceptAll = zeros(nChains)
        ProgBar = Progress(nChains*(chunkSize-1), dt=1, desc="Running chunk $ii", barlen=50, enabled=verbose)
        Threads.@threads for jj = 1:nChains
            initialValues = missing
            if ii > 1
                initialValues = θall[firstInd-1,:,jj]
            end
            (θall[inds,:,jj], nAcceptAll[jj]) = runSingleEpochPMMH(bootstrapFilter, covEst, Y, opts; ProgBar=ProgBar, chainNum=jj, initialValues=initialValues)
        end
        finish!(ProgBar)
        
        # While in the tuning phase....
        if isTuning
            
            # Get chains object
            chains = Chains(θall[inds,:,:])
            dfc = DataFrame(summarize(chains))
            
            # Upate mean est and covariance matrix
            frobeniusOld = norm(covEst)
            θall_tmp = resizeParams(θall[inds,:,:])
            (meanEst, covEst) = updateCovariance(meanEst, covEst, θall_tmp, (100*(ii-1))+1, chunkSize) # Uses n_old = 100 in the first iteration to prevent stochasticity overruling
            frobeniusNew = norm(covEst)
            change = frobeniusNew / frobeniusOld
            if verbose
                println("Old frobenius: $frobeniusOld, new frobenius: $frobeniusNew, change: $change")
                println("Approx covariance matrix of parameter samples: $covEst")
            end

            # Save diagnostics
            if maximum(change) < 1.25 && minimum(change) > 0.8
                isTuning = false
                maxTuningInd = lastInd
            end
            
            diag_tmp = DataFrame()
            diag_tmp.chunk = repeat([ii], nParams)
            diag_tmp.tuning = repeat([true], nParams)
            diag_tmp.pAccept = repeat([sum(nAcceptAll)/(nChains * (chunkSize - 1))],  nParams)
            diag_tmp.Rhat = dfc.rhat
            diag_tmp.ESS = dfc.ess
            diag_tmp.mean = dfc.mean
            diag_tmp.std = dfc.std
            diag = vcat(diag, diag_tmp)
            
            if verbose
                println(diag_tmp)
            end

            nTuningChunks = nTuningChunks + 1
            
        else
            
            st_ind = Int(round(maxTuningInd + ((ii - nTuningChunks) * chunkSize * 0.2))) # Ignore 20% of a chunk from the start each time
            chains = Chains(θall[st_ind:lastInd,:,:])
            dfc = DataFrame(summarize(chains))
            
            # Update mean est and covariance matrix
            θall_tmp = resizeParams(θall[inds,:,:])
            (meanEst, covEst) = updateCovariance(meanEst, covEst, θall_tmp, 100*ii, chunkSize)
            if verbose
                println("Frobenius of covariance matrix: $(norm(covEst))")
                println("Approx covariance matrix of parameter samples: $covEst")
            end

            # Fetch max Rhat and min ESS
            Rhat = maximum(dfc.rhat)
            ESS = minimum(dfc.ess)
            
            # Save diagnostics
            diag_tmp = DataFrame()
            diag_tmp.chunk = repeat([ii], nParams)
            diag_tmp.tuning = repeat([false], nParams)
            diag_tmp.pAccept = repeat([sum(nAcceptAll)/(nChains * (chunkSize - 1))],  nParams)
            diag_tmp.Rhat = dfc.rhat
            diag_tmp.ESS = dfc.ess
            diag_tmp.mean = dfc.mean
            diag_tmp.std = dfc.std
            diag = vcat(diag, diag_tmp)
            
            if verbose
                println(diag_tmp)
            end
            
        end
        
        # Increment the number of chunks
        ii = ii + 1
        
    end
    
    if isTuning
        @warn("Algorithm never finished tuning. Returning PMMH results with all samples.")
        maxTuningInd = 0
    end
    
    st_ind = Int(round(maxTuningInd + (((ii - 1) - nTuningChunks) * chunkSize * 0.2)))
    lastInd = (ii-1)*chunkSize
    return(θall[st_ind:lastInd, :, :], diag)
    
end








function updateCovariance(μ_old, Σ_old, X_batch, n_old, n_batch)

    # Calculate batch mean and covariance
    μ_batch = mean(X_batch, dims=1)
    Σ_batch = cov(X_batch)

    # Update mean
    μ_new = (n_old * μ_old + n_batch * μ_batch) / (n_old + n_batch)

    # Update covariance
    mean_diff = μ_batch - μ_old
    Σ_new = (n_old * Σ_old + n_batch * Σ_batch + (n_old * n_batch) / (n_old + n_batch) * (mean_diff' * mean_diff)) / (n_old + n_batch)

    return μ_new, Σ_new
end



function runSingleEpochPMMH(bootstrapFilter::Function, covEst, Y::DataFrame, opts::Dict; ProgBar=missing, chainNum=0, initialValues=missing, verbose=false)
    
    # Extract options that are frequently used
    chunkSize = opts["chunkSize"]
    nParams = length(opts["paramPriors"])
    
    # Pre-allocate outputs
    θ = zeros(chunkSize, nParams)
    loglik = -Inf
    nAccept = 0
    
    # Sample initial parameter values and compute initial log-likelihood (only if no previous values are included)
    if ismissing(initialValues)
        (θ[1,:], loglik) = sampleInitialValues(bootstrapFilter, Y, opts)
    else
        θ[1,:] = initialValues
        loglik = estimateLogLik(bootstrapFilter, initialValues, Y, opts; ignoreerror=true)
    end

    # Setup progress bar if one is not provided
    if ismissing(ProgBar)
        ProgBar = Progress(chunkSize-1, dt=1, desc="Running single chunk...", barlen=50, enabled=verbose)
    end
        
    # Run PMMH
    for ii = 2:chunkSize

        # Propose new parameter value
        θprop = proposeParameterValue(θ[ii-1,:], nParams, covEst, opts)
        
        # Estimate log-likelihood of the proposal
        loglikProp = estimateLogLik(bootstrapFilter, θprop, Y, opts; ignoreerror=true)
        
        # Calculate acceptance probability
        α = calculateAcceptanceProbability(loglikProp, loglik, θprop, θ[ii-1,:], covEst, opts)
        
        # Accept or reject proposal
        if rand() < α
            θ[ii,:] = θprop
            loglik = loglikProp
            nAccept += 1
        else
            θ[ii,:] = θ[ii-1,:]
        end
        
        next!(ProgBar)
        
    end
    
    return(θ, nAccept)
    
end



function sampleInitialValues(bootstrapFilter::Function, Y::DataFrame, opts::Dict)

    initalValuesSampled = false
    nAttempts = 0
    θ = missing
    loglik = missing

    while !initalValuesSampled
        nAttempts = nAttempts + 1
        θ = [rand(opts["initialParamSamplers"][ii]) for ii in 1:length(opts["initialParamSamplers"])]
        try
            loglik = estimateLogLik(bootstrapFilter, θ, Y, opts)
            initalValuesSampled = true
        catch e
            println("Initial value θ = $θ failed. Trying again...")
        end
        if nAttempts > 100
            error("Error sampling initial values. Either initial values could not be sampled, or there is a bug in the code.")
        end
    end

    return(θ, loglik)

end


function proposeParameterValue(θold, nParams, covEst, opts; maxAttempts=100000)

    counter = 1
    validParamFound = true
    violatingInd = missing

    # Sample first attempt and check
    θprop = rand(MvNormal(θold, ((2.38^2)/nParams) * covEst))
    for kk in 1:nParams
        if θprop[kk] < opts["paramLimits"][kk][1] || θprop[kk] > opts["paramLimits"][kk][2]
            validParamFound = false
            violatingInd = kk
        end
    end

    # If the first attempt is invalid, keep sampling until a valid parameter is found
    while !validParamFound & (counter <= maxAttempts)
        validParamFound = true
        θprop = rand(MvNormal(θold, ((2.38^2)/nParams) * covEst))
        for kk in 1:nParams
            if θprop[kk] < opts["paramLimits"][kk][1] || θprop[kk] > opts["paramLimits"][kk][2]
                validParamFound = false
                violatingInd = kk
            end
        end
        counter += 1
    end

    # If we hit the maxAttemps, throw an error
    if counter >= maxAttempts
        error("Could not sample a valid parameter value starting from $θold.  Final attempt was $θprop at index $violatingInd. Check the parameter limits.")
    end

    return(θprop)

end


function calculateAcceptanceProbability(loglikprop, loglikprev, θprop, θprev, covEst, opts)
    
    nParams = length(θprop)
    
    # Likelihood ratio
    likelihoodratio = exp(loglikprop - loglikprev)
    
    priorratio = missing
    try
        priorratio = prod([pdf(opts["paramPriors"][jj], θprop[jj]) / pdf(opts["paramPriors"][jj], θprev[jj]) for jj in 1:nParams])
    catch e
        println("ERROR IN PRIOR RATIO CALCULATION")
        println("θprop: $θprop")
        println("θprev: $θprev")
        rethrow(e)
    end

    proposalratio = pdf(MvNormal(θprev, covEst), θprop) / pdf(MvNormal(θprop, covEst), θprev)
    
    acceptanceprob = min(1, likelihoodratio * priorratio * proposalratio)
    
    return(acceptanceprob)
    
end
