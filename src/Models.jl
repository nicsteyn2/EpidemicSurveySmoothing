# We redefine some PDF functions here for speed, these are only used in more complicated models such as the variants model
include("Support.jl")

# This file contains the SMC models used throughout the project:
#
# Main paper:
#  - BinomialModel()
#  - BetaBinomialModel()
#  - WeightedModel()
#  - RtModel() (and the ever so slightly faster RtModelForPaper())
#
# Supplementary material:
#  - LogitModel()
#  - VariantsModel8to13()
#  - VariantsModel14to19()


# Basic model (binomial observation dist)
function BinomialModel(θ::Vector, Y::DataFrame, opts::Dict; showProgress=false)
    
    T = nrow(Y) # Length of data
    N = opts["N"] # Number of particles
    
    # Pre-allocate outputs
    R = zeros(N, T) # Hidden state for growth rate r_t
    P = zeros(N, T) # Hidden state for true prevalence P_t
    W = zeros(N, T) # Weights
    
    # Initialise states
    R[:,1] = rand(opts["pr0"], N)
    P[:,1] = rand(opts["pP0"], N)
    
    # Run the filter
    ProgBar = Progress(T-1, dt=1, desc="Running particle filter...", barlen=50, enabled=showProgress)
    for tt in 2:T
        
        # Project hidden-states
        R[:,tt] = rand.(Normal.(R[:,tt-1], θ[1]))
        P[:,tt] = clamp.(P[:,tt-1] .* exp.(R[:,tt]), 0, 1)
        
        # Weight and resample, but only if there are observed data
        if Y.nSamples[tt] > 0
            # Weight
            W[:,tt] = pdf.(Binomial.(Y.nSamples[tt], P[:,tt]), Y.nPos[tt])
            # Resample
            inds = wsample(1:N, W[:,tt], N)
            st_time = max(tt - opts["L"], 1) # Fixed-lag resampling
            R[:,st_time:tt] = R[inds,st_time:tt]
            P[:,st_time:tt] = P[inds,st_time:tt]
        else
            W[:,tt] .= 1   
        end
        
        # Update progress bar
        next!(ProgBar)
        
    end
    
    # Store output in a consistent format
    X = zeros(N, T, 2)
    X[:,:,1] = R
    X[:,:,2] = P
    
    return(X, W)
    
end



# Extra-binomial variation (beta-binomial observation dist)
function BetaBinomialModel(θ::Vector, Y::DataFrame, opts::Dict; showProgress=false)
    
    T = nrow(Y) # Length of data
    N = opts["N"] # Number of particles
    
    # Pre-allocate outputs
    R = zeros(N, T) # Hidden state for growth rate r_t
    P = zeros(N, T) # Hidden state for true prevalence P_t
    W = zeros(N, T) # Weights
    
    # Initialise states
    R[:,1] = rand(opts["pr0"], N)
    P[:,1] = rand(opts["pP0"], N)
    
    # Run the filter
    ProgBar = Progress(T-1, dt=1, desc="Running particle filter...", barlen=50, enabled=showProgress)
    for tt in 2:T
        
        # Project hidden-states
        R[:,tt] = rand.(Normal.(R[:,tt-1], θ[1]))
        P[:,tt] = clamp.(P[:,tt-1] .* exp.(R[:,tt]), 0, 1)
        
        # Weight and resample, but only if there are observed data
        if Y.nSamples[tt] > 0
            # Weight
            α = P[:,tt] .* ( (1/θ[2]) - 1)
            β = (1 .- P[:,tt]) * ( (1/θ[2]) - 1)
            W[:,tt] = pdf.(BetaBinomial.(Y.nSamples[tt], α, β), Y.nPos[tt])
            # Resample
            inds = wsample(1:N, W[:,tt], N)
            st_time = max(tt - opts["L"], 1) # Fixed-lag resampling
            R[:,st_time:tt] = R[inds,st_time:tt]
            P[:,st_time:tt] = P[inds,st_time:tt]
        else
            W[:,tt] .= 1   
        end
        
        # Update progress bar
        next!(ProgBar)
        
    end
    
    # Store output in a consistent format
    X = zeros(N, T, 2)
    X[:,:,1] = R
    X[:,:,2] = P
    
    return(X, W)
    
end



# Weighted model
function WeightedModel(θ::Vector, Y::DataFrame, opts::Dict; showProgress=false)
    
    T = nrow(Y) # Length of data
    N = opts["N"] # Number of particles
    
    # Pre-allocate outputs
    R = zeros(N, T) # Hidden state for growth rate r_t
    P = zeros(N, T) # Hidden state for true prevalence P_t
    W = zeros(N, T) # Weights
    
    # Initialise states
    R[:,1] = rand(opts["pr0"], N)
    P[:,1] = rand(opts["pP0"], N)
    
    # Run the filter
    ProgBar = Progress(T-1, dt=1, desc="Running particle filter...", barlen=50, enabled=showProgress)
    for tt in 2:T
        
        # Calculate rt bounds if used (useful when true prevalence and nt are both low)
        rt_min = -Inf
        rt_max = Inf
        if opts["boundPrevalence"]
            rt_min = log(opts["minPrevalence"]) .- log.(P[:,tt-1])
            rt_max = log(opts["maxPrevalence"]) .- log.(P[:,tt-1])
        end
        
        # Project hidden-states
        R[:,tt] = rand.(Truncated.(Normal.(R[:,tt-1], θ[1]), rt_min, rt_max))
        P[:,tt] = clamp.(P[:,tt-1] .* exp.(R[:,tt]), 0, 1)
        
        # Weight and resample, but only if there are observed data
        if Y.nSamples[tt] > 0
            
            # Weight
            phat_var = θ[2] * Y.ssw[tt] * P[:,tt] .* (1 .- P[:,tt])
            W[:,tt] = pdf.(Normal.(P[:,tt], sqrt.(phat_var)), Y.phat[tt])
            
            # Resample hidden-states
            inds = wsample(1:N, W[:,tt], N)
            st_time = max(tt - opts["L"], 1)
            R[:,st_time:tt] = R[inds,st_time:tt]
            P[:,st_time:tt] = P[inds,st_time:tt]
        else
            W[:,tt] .= 1   
        end
        
        # Update progress bar
        next!(ProgBar)
        
    end
    
    # Store output in a consistent format
    X = zeros(N, T, 2)
    X[:,:,1] = R
    X[:,:,2] = P
    
    return(X, W)
    
end