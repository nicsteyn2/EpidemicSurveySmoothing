
using Random, Distributions, ProgressMeter, DataFrames

# ---------------------------------------------------------------------
# ---------------- Epidemic simulator functions -----------------------
# ---------------------------------------------------------------------

function simulateRandomWalk( ; σ = 0.01, P0 = 0.01, r0 = 0, tMax = 100, rt_max=0.3, rt_min=-0.3)
    
    # Pre-allocate output
    rt = zeros(tMax)
    Pt = zeros(tMax)
    
    # Simulate the epidemic
    tt = 1
    rt[tt] = rand(Truncated(Normal(r0, σ), rt_min, rt_max))
    Pt[tt] = P0 * exp(rt[tt])
    
    for tt = 2:tMax
        rt[tt] = rand(Truncated(Normal(rt[tt-1], σ), rt_min, rt_max))
        Pt[tt] = Pt[tt-1] * exp(rt[tt])
    end
    
    return(rt, Pt)
    
end

function simulateBoundedRandomWalk(; Pt_min=0.001, Pt_max=0.3, maxTrials=1e5, σ=0.01, P0=0.01, r0=0, tMax=100, rt_max=0.3, rt_min=-0.3)

    # Simulate first epidemic
    (rt, Pt) = simulateRandomWalk(σ=σ, P0=P0, r0=r0, tMax=tMax, rt_max=rt_max, rt_min=rt_min)
    counter = 1

    # Check bounds
    finished = (maximum(Pt) < Pt_max) & (minimum(Pt) > Pt_min)

    # Loop until we find a solution or reach the maximum number of trials
    while (counter <= maxTrials && !finished)
        (rt, Pt) = simulateRandomWalk(σ=σ, P0=P0, r0=r0, tMax=tMax, rt_max=rt_max, rt_min=rt_min)
        counter += 1
        finished = (maximum(Pt) < Pt_max) & (minimum(Pt) > Pt_min)
    end

    # Return the solution if we found one
    if finished
        return(rt, Pt)
    else
        error("Failed to simulate an epidemic within the specified bounds.")
    end

end



# #TODO: Update with the new sampling regime
# function simulateAR1(; ϕ = 0.9, σ = 0.01, P0 = 0.01, r0 = 0, tMax = 100, Pt_max=0.3, Pt_min=0.005, rt_max=0.3, rt_min=-0.3)

#     rt = zeros(tMax)
#     Pt = zeros(tMax)

#     # Simulate the epidemic
#     tt = 1
#     rt[tt] = rand(Truncated(Normal(r0, σ), rt_min, rt_max))
#     Pt[tt] = P0 * exp(rt[tt])

#     for tt = 2:tMax
#         rt_max_t = min(log(Pt_max/Pt[tt-1]), rt_max)
#         rt_min_t = max(log(Pt_min/Pt[tt-1]), rt_min)
#         rt[tt] = rand(Truncated(Normal(ϕ * rt[tt-1], σ), rt_min_t, rt_max_t))
#         Pt[tt] = Pt[tt-1] * exp(rt[tt])
#     end

#     return(rt, Pt)

# end



# ---------------------------------------------------------------------
# ------------------ Survey simulation functions ----------------------
# ---------------------------------------------------------------------

function simulateSurvey(nt, rt, Pt; method="binomial", params=[])

    Y = DataFrame()
    Y.obsmethod = repeat([method], length(Pt))
    Y.t = 1:length(Pt)
    Y.rt = rt
    Y.Pt = Pt
    Y.nSamples = repeat([nt], length(Pt))

    if method == "binomial"
        Y.nPos = rand.(Binomial.(Y.nSamples, Pt))
        Y.phat = Y.nPos ./ Y.nSamples
        Y.ssw = 1 ./ Y.nSamples
    elseif method == "betabinom"
        ρ = params[1]
        if ρ == 0
            Y.nPos = rand.(Binomial.(Y.nSamples, Pt))
        else
            αt = Pt .* (1/ρ - 1)
            βt = (1 .- Pt) .* (1/ρ - 1)
            Y.nPos = rand.(BetaBinomial.(Y.nSamples, αt, βt))
        end
        Y.phat = Y.nPos ./ Y.nSamples
        Y.ssw = 1 ./ Y.nSamples
    elseif method == "weighted"
        ξ = params[1]
        if ξ == 0
            Y.nPos = rand.(Binomial.(Y.nSamples, Pt))
            Y.phat = Y.nPos ./ Y.nSamples
            Y.ssw = 1 ./ Y.nSamples
            Y.theoreticalvar = repeat([NaN], length(Pt))
        else
            Y.nPos = zeros(Int, length(Pt))
            Y.phat = zeros(length(Pt))
            Y.ssw = zeros(length(Pt))
            Y.theoreticalvar = zeros(length(Pt))
            for tt = 1:length(Pt)
                Wraw = rand(LogNormal(0, ξ), Y.nSamples[tt])
                W = Wraw/sum(Wraw)
                Pi = clamp.(Pt[tt] * (W/sum(W .^ 2)), 0, 1)
                Xi = rand.(Bernoulli.(Pi))
                Y.nPos[tt] = sum(Xi)
                Y.phat[tt] = sum(W .* Xi)
                Y.ssw[tt] = sum(W .^ 2)
                Y.theoreticalvar[tt] = sum( (W .^ 2) .* Pi .* (1 .- Pi) )
            end
        end
    else
        error("Survey simulation method not recognised")
    end

    return(Y)

end