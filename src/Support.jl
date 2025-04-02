
using Distributions



function processResults(X)

    m = vec(mean(X, dims=1))
    med = [quantile(Xt, 0.5) for Xt in eachcol(X)]
    l = [quantile(Xt, 0.025) for Xt in eachcol(X)]
    u = [quantile(Xt, 0.975) for Xt in eachcol(X)]

    return(m, med, l, u)

end


function resizeParams(θ)
    
    if length(size(θ)) != 3
        error()
    end
    
    nSamples = size(θ)[1]
    nParams = size(θ)[2]
    nChains = size(θ)[3]
    
    θout = zeros(nSamples * nChains, nParams)
    for ii = 1:nChains
        inds = ((ii-1)*nSamples + 1):(ii*nSamples)
        θout[inds,:] = θ[:,:,ii]
    end
    
    return(θout)
    
end


# Fast distributions

import SpecialFunctions: loggamma
import LogExpFunctions: xlogy

function fastNegativeBinomialPDF(x::Real, r::AbstractVector, p::Real)

    constant_term = -loggamma(x + 1) + xlogy(x, 1-p)
    variable_term = loggamma.(x .+ r) - loggamma.(r) + xlogy.(r, p)
    if x > 0
        return exp.(constant_term .+ variable_term)
    else
        return ifelse.(r.==0, 1, exp.(constant_term .+ variable_term))
    end

end


function fastMultinomialLogPDF(n::Int, x::Vector{Int}, p::Matrix)
    
    common_values = loggamma(n+1) - sum(loggamma.(x .+ 1))
    individual_values = [sum(xlogy.(x, p[ii,:])) for ii = 1:size(p,1)]
    return individual_values .+ common_values
    
end