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



function RtModel(θ::Vector, Y::DataFrame, opts::Dict; showProgress=false)
    
    # Extract options that are frequently used
    T = nrow(Y)
    N = opts["N"]
    init = opts["initPeriod"]
    
    # Extract fixed parameters and distributions
    genTimeDist = opts["genTimeDist"]
    infToPos = opts["infToPositivity"]
    
    # Pre-allocate outputs
    R = zeros(N, T+init)
    I = zeros(N, T+init)
    P = zeros(N, T+init)
    W = zeros(N, T+init)
    
    # Initialise states
    R[:,init] .= rand(opts["pR0"], N)
    P[:,1:init] .= rand(opts["pP0"], N)
    I[:,1:init] = P[:,1:init] / sum(infToPos)
    
    # Run the filter
    ProgBar = Progress(T-1, dt=1, desc="Running particle filter...", barlen=50, enabled=showProgress)
    for tt in (init+1):(T+init)
        
        tNow = tt - init # tt is the index of the hidden-state matrices, tNow is the corresponding index of the data
        
        # Calculate force-of-infection
        st_ind = max(1, tt - opts["L"] + 1)
        en_ind = min(tt, opts["L"])
        FOI = sum(I[:, tt:-1:st_ind] .* genTimeDist[1:1:en_ind]', dims=2)/sum(genTimeDist[1:1:en_ind])
        
        # Project hidden-states
        R[:,tt] = exp.(rand.(Normal.(log.(R[:,tt-1]), θ[1])))
        I[:,tt] = R[:,tt] .* FOI
        scalingFactor = sum(infToPos)/sum(infToPos[1:1:en_ind]) # Account for init periods shorter than observation times
        P[:,tt] = clamp.(sum(scalingFactor * I[:,tt:-1:st_ind] .* infToPos[1:1:en_ind]', dims=2), 0, 1)
        
        # Weight and resample
        if Y.nSamples[tNow] > 0
            
            # Weight hidden states
            α = P[:,tt] .* ( (1/θ[2]) - 1)
            β = (1 .- P[:,tt]) * ( (1/θ[2]) - 1)
            W[:,tt] = pdf.(BetaBinomial.(Y.nSamples[tNow], α, β), Y.nPos[tNow])
            
            # Resample
            inds = wsample(1:N, W[:,tt], N)
            st_time = max(tt - opts["L"], 1) # Fixed-lag resampling
            R[:,st_time:tt] = R[inds,st_time:tt]
            I[:,st_time:tt] = I[inds,st_time:tt]
            P[:,st_time:tt] = P[inds,st_time:tt]
            
        else
            
            W[:,tt] .= 1
            
        end
        
        next!(ProgBar)
        
    end
    
    # Store output in a consistent format
    X = zeros(N, T, 3)
    X[:,:,1] = R[:,init+1:end]
    X[:,:,2] = I[:,init+1:end]
    X[:,:,3] = P[:,init+1:end]
    
    return(X, W[:,init+1:end])
    
end


function RtModelForPaper(θ::Vector, Y::DataFrame, opts::Dict; showProgress=false)

    # A slightly modified version that assumes the function is only being used for PMMH if L = 30
    # This allows it to skip some of the resampling steps, increasing speed a little
    
    # Extract options that are frequently used
    T = nrow(Y)
    N = opts["N"]
    init = opts["initPeriod"]
    
    # Extract fixed parameters and distributions
    genTimeDist = opts["genTimeDist"]
    infToPos = opts["infToPositivity"]
    
    # Pre-allocate outputs
    R = zeros(N, T+init)
    I = zeros(N, T+init)
    P = zeros(N, T+init)
    W = zeros(N, T+init)
    
    # Initialise states
    R[:,init] .= rand(opts["pR0"], N)
    P[:,1:init] .= rand(opts["pP0"], N)
    I[:,1:init] = P[:,1:init] / sum(infToPos)
    
    # Run the filter
    ProgBar = Progress(T-1, dt=1, desc="Running particle filter...", barlen=50, enabled=showProgress)
    for tt in (init+1):(T+init)
        
        tNow = tt - init # tt is the index of the hidden-state matrices, tNow is the corresponding index of the data
        
        # Calculate force-of-infection
        st_ind = max(1, tt - opts["L"] + 1)
        en_ind = min(tt, opts["L"])
        FOI = sum(I[:, tt:-1:st_ind] .* genTimeDist[1:1:en_ind]', dims=2)/sum(genTimeDist[1:1:en_ind])
        
        # Project hidden-states
        R[:,tt] = exp.(rand.(Normal.(log.(R[:,tt-1]), θ[1])))
        I[:,tt] = R[:,tt] .* FOI
        scalingFactor = sum(infToPos)/sum(infToPos[1:1:en_ind]) # Account for init periods shorter than observation times
        P[:,tt] = clamp.(sum(scalingFactor * I[:,tt:-1:st_ind] .* infToPos[1:1:en_ind]', dims=2), 0, 1)
        
        # Weight and resample
        if Y.nSamples[tNow] > 0
            
            # Weight hidden states
            α = P[:,tt] .* ( (1/θ[2]) - 1)
            β = (1 .- P[:,tt]) * ( (1/θ[2]) - 1)
            W[:,tt] = pdf.(BetaBinomial.(Y.nSamples[tNow], α, β), Y.nPos[tNow])
            
            # Resample
            inds = wsample(1:N, W[:,tt], N)
            st_time = max(tt - opts["L"], 1) # Fixed-lag resampling
            if opts["L"] == 30 # We only use L = 30 in the paper inside PMMH, in which case we can skip much of the resampling
                R[:,tt] = R[inds,tt]
                I[:,st_time:tt] = I[inds,st_time:tt]
            else
                R[:,st_time:tt] = R[inds,st_time:tt]
                I[:,st_time:tt] = I[inds,st_time:tt]
                P[:,st_time:tt] = P[inds,st_time:tt]
            end
            
        else
            
            W[:,tt] .= 1
            
        end
        
        next!(ProgBar)
        
    end
    
    # Store output in a consistent format
    X = zeros(N, T, 3)
    X[:,:,1] = R[:,init+1:end]
    X[:,:,2] = I[:,init+1:end]
    X[:,:,3] = P[:,init+1:end]
    
    return(X, W[:,init+1:end])
    
end