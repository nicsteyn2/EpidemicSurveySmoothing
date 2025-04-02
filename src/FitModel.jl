
using DataFrames

include("Support.jl")
include("PMMH.jl")
include("MarginalPosterior.jl")
include("Likelihood.jl")

# This function performs the start-to-finish model fitting procedure
function fitModel(bootstrapFilter::Function, Y::DataFrame, optsIn::Dict; skipResamplingPMMH=false, verbose=false, posteriorPredictive="none", returnSamples=false, checkLogLikStdDev=false)

    # Deepcopy the options to prevent unexpected modification
    opts = deepcopy(optsIn)
    
    # Extract options
    paramNames = opts["paramNames"] # PMMH parameter names (useful for storing tidy output)
    stateNames = opts["stateNames"] # Hidden-state names
    nParams = length(opts["paramPriors"])
    θtest = missing
    if haskey(opts, "initialParamSamplers")
        θtest = [mean(opts["initialParamSamplers"][ii]) for ii in 1:nParams]
    else
        θtest = [mean(opts["paramPriors"][ii]) for ii in 1:nParams]
    end
    
    # If the model is Markovian, we can skip resampling when estimating the likelihood
    L = opts["L"] # Store this to use when finding the marginal posterior distribution
    if skipResamplingPMMH
        opts["L"] = 1
    end
    
    # Check the standard deviation of the log-likelihood
    if checkLogLikStdDev
        (loglikstddev, _) = estimateStdDevLogLik(100, bootstrapFilter, θtest, Y, opts; showProgress=verbose)
        if verbose
            println("Estimated standard deviation of the log-likelihood: ", loglikstddev)
        end
        if loglikstddev > 3
            error("The standard deviation of the log-likelihood is too high. Consider increasing the number of particles used, or changing the model.")
        elseif loglikstddev > 2
            @warn("The standard deviation of the log-likelihood is high (>2). Consider increasing the number of particles used, or changing the model.")
        end
    end
    
    # Start the timer
    S = time()
    # Run PMMH
    (θ, diag) = PMMH(bootstrapFilter, Y, opts; verbose=verbose)
    chains = Chains(θ, opts["paramNames"])
    
    # Store parameter values
    df_params = DataFrame()
    df_params.param = paramNames
    df_params.m = [mean(θ[:,ii,:]) for ii in 1:size(θ)[2]]
    df_params.l = [quantile(vec(θ[:,ii,:]), 0.025) for ii in 1:size(θ)[2]]
    df_params.u = [quantile(vec(θ[:,ii,:]), 0.975) for ii in 1:size(θ)[2]]
    df_params.rhat = DataFrame(summarize(chains)).rhat
    df_params.ess = DataFrame(summarize(chains)).ess
    
    # Reset resampling window
    if "L-marginalposterior" in keys(opts)
        opts["L"] = opts["L-marginalposterior"]
    else
        opts["L"] = L
    end
    
    # Fetch samples from the marginal posterior distribution
    (X, θs) = marginalPosterior(bootstrapFilter, θ, Y, opts)
    
    # Extract dimensions
    nStates = missing
    if length(size(X)) == 2
        nStates = 1
    else
        nStates = size(X)[3]
    end
    
    # Process marginal posterior samples
    df_states = DataFrame()
    if nStates == 1
        (m, _, l, u) = processResults(X)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable=stateNames[1])
        df_states = vcat(df_states, tmp)
    else
        for ii = 1:nStates
            stateName = stateNames[ii]
            (m, _, l, u) = processResults(X[:,:,ii])
            tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable=stateName)
            df_states = vcat(df_states, tmp)
        end
    end
    
    # Fetch posterior predictive distribution
    if posteriorPredictive == "binomial"
        P = X[:,:,end] # Always assume the last state is prevalence
        Xpred = similar(P)
        for tt = 1:nrow(Y)
            Xpred[:,tt] = rand.(Binomial.(Y.nSamples[tt], P[:,tt]))
        end
        (m, _, l, u) = processResults(Xpred)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable="nPos")
        tmp2 = DataFrame(t=Y.t, mean=m ./ Y.nSamples, lower=l ./ Y.nSamples, upper=u ./ Y.nSamples, variable="Obs positivity")
        df_states = vcat(df_states, tmp)
        df_states = vcat(df_states, tmp2)
    elseif posteriorPredictive == "betabinom"
        P = X[:,:,end] # Always assume the last state is prevalence
        Xpred = similar(P)
        for tt = 1:nrow(Y)
            α = P[:,tt] .* ( (1 ./ θs[:,2]) .- 1)
            β = (1 .- P[:,tt]) .* ( (1 ./ θs[:,2]) .- 1)
            Xpred[:,tt] = rand.(BetaBinomial.(Y.nSamples[tt], α, β))
        end
        (m, _, l, u) = processResults(Xpred)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable="nPos")
        tmp2 = DataFrame(t=Y.t, mean=m ./ Y.nSamples, lower=l ./ Y.nSamples, upper=u ./ Y.nSamples, variable="Obs positivity")
        df_states = vcat(df_states, tmp)
        df_states = vcat(df_states, tmp2)
    elseif posteriorPredictive == "weighted"
        P = X[:,:,end]
        Xpred = similar(P)
        for tt = 1:nrow(Y)
            phat_var = θs[:,2] .* Y.ssw[tt] .* P[:,tt] .* (1 .- P[:,tt])
            Xpred[:,tt] = rand.(Normal.(P[:,tt], sqrt.(phat_var)))
        end
        (m, _, l, u) = processResults(Xpred)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable="phat")
        df_states = vcat(df_states, tmp)
    end
    
    # If Y contains a date vector, include it in our outputs
    if "date" in names(Y)
        df_states = leftjoin(df_states, Y[:,[:t, :date]], on=:t)
    end
    
    # End the timer
    E = time()
    if verbose
        println("Model finished. Time elapsed: ", round(E-S, sigdigits=3), " seconds")
    end
    
    if returnSamples && posteriorPredictive == "none"
        return(df_states, df_params, θ, diag, X)
    elseif returnSamples
        return(df_states, df_params, θ, diag, X, Xpred)
    else
        return(df_states, df_params, θ, diag)
    end

end


# This function performs the start-to-finish model fitting procedure, but with fixed parameters (i.e. no PMMH)
function fitModelFixedParams(θ::Vector, bootstrapFilter::Function, Y::DataFrame, opts::Dict; posteriorPredictive="none")


    # Run bootstrap filter
    (X, _) = bootstrapFilter(θ, Y, opts)
    if length(size(X)) == 2
        nStates = 1
    else
        nStates = size(X)[3]
    end

    # Process samples
    df_states = DataFrame()
    if nStates == 1
        (m, _, l, u) = processResults(X)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable=opts["stateNames"][1])
        df_states = vcat(df_states, tmp)
    else
        for ii = 1:nStates
            stateName = opts["stateNames"][ii]
            (m, _, l, u) = processResults(X[:,:,ii])
            tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable=stateName)
            df_states = vcat(df_states, tmp)
        end
    end

    # Fetch posterior predictive distribution
    if posteriorPredictive == "binomial"
        P = X[:,:,2]
        Xpred = similar(P)
        for tt = 1:nrow(Y)
            Xpred[:,tt] = rand.(Binomial.(Y.nSamples[tt], P[:,tt]))
        end
        (m, _, l, u) = processResults(Xpred)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable="nPos")
        tmp2 = DataFrame(t=Y.t, mean=m ./ Y.nSamples, lower=l ./ Y.nSamples, upper=u ./ Y.nSamples, variable="Obs positivity")
        df_states = vcat(df_states, tmp)
        df_states = vcat(df_states, tmp2)
    elseif posteriorPredictive == "betabinom"
        P = X[:,:,2]
        Xpred = similar(P)
        for tt = 1:nrow(Y)
            α = P[:,tt] .* ( (1 ./ θ[2]) .- 1)
            β = (1 .- P[:,tt]) .* ( (1 ./ θ[2]) .- 1)
            Xpred[:,tt] = rand.(BetaBinomial.(Y.nSamples[tt], α, β))
        end
        (m, _, l, u) = processResults(Xpred)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable="nPos")
        tmp2 = DataFrame(t=Y.t, mean=m ./ Y.nSamples, lower=l ./ Y.nSamples, upper=u ./ Y.nSamples, variable="Obs positivity")
        df_states = vcat(df_states, tmp)
        df_states = vcat(df_states, tmp2)
    elseif posteriorPredictive == "weighted"
        P = X[:,:,2]
        Xpred = similar(P)
        for tt = 1:nrow(Y)
            phat_var = θ[2] .* Y.ssw[tt] .* P[:,tt] .* (1 .- P[:,tt])
            Xpred[:,tt] = rand.(Normal.(P[:,tt], sqrt.(phat_var)))
        end
        (m, _, l, u) = processResults(Xpred)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable="Obs positivity")
        df_states = vcat(df_states, tmp)
    end
    
    # If Y contains a date vector, include it in our outputs
    if "date" in names(Y)
        df_states = leftjoin(df_states, Y[:,[:t, :date]], on=:t)
    end

    return(df_states)

end
