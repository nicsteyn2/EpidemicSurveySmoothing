
using CSV, DataFrames


include("../src/Models.jl")
include("../src/FitModel.jl")
include("../src/Opts.jl")


# # Simulate and store a simulation (if the data does not already exist)
# include("../src/Simulations.jl")
# (rt, Pt) = simulateBoundedRandomWalk()
# Y = simulateSurvey(5000, rt, Pt; method="binomial")
# CSV.write("paper/outputs/s1-sensitivity/s1_exampledata.csv", Y)


# Load the example data (if it exists)
Y = CSV.read("paper/outputs/s1-sensitivity/s1_exampledata.csv", DataFrame)



# ------------ TEST INITIAL DISTRIBUTION OF r0 ------------

# Fetch opts
opts = fetchDefaultOpts("binomial")

# Specify the initial distribution for Pt
phat = max(Y.nPos[1], 1)/Y.nSamples[1]
sd = sqrt((phat*(1-phat))/Y.nSamples[1])
opts["pP0"] = Truncated(Normal(phat, sd), 0, 1)

df_states_all = DataFrame()
df_params_all = DataFrame()

opts["pr0"] = Uniform(-0.1, 0.1)
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Uniform(-0.1, 0.1) (narrower)"
df_params[!,:prior] .= "Uniform(-0.1, 0.1) (narrower)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

opts["pr0"] = Uniform(-0.3, 0.3)
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Uniform(-0.3, 0.3) (default)"
df_params[!,:prior] .= "Uniform(-0.3, 0.3) (default)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

opts["pr0"] = Uniform(-1, 1)
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Uniform(-1.0, 1.0) (wider)"
df_params[!,:prior] .= "Uniform(-1.0, 1.0) (wider)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

CSV.write("paper/outputs/s1-sensitivity/s1_sensitivity_pr0_states.csv", df_states_all)
CSV.write("paper/outputs/s1-sensitivity/s1_sensitivity_pr0_params.csv", df_params_all)


# ------------ TEST INITIAL DISTRIBUTION OF P0 ------------

# Fetch opts
opts = fetchDefaultOpts("binomial")

# Specify the initial distribution for Pt
phat = max(Y.nPos[1], 1)/Y.nSamples[1]
sd = sqrt((phat*(1-phat))/Y.nSamples[1])
opts["pP0"] = Truncated(Normal(phat, sd), 0, 1)

df_states_all = DataFrame()
df_params_all = DataFrame()

opts["pP0"] = Truncated(Normal(phat, sd), 0, 1)
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Wald interval (default)"
df_params[!,:prior] .= "Wald interval (default)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

opts["pP0"] = Truncated(Normal(phat, 2*sd), 0, 1)
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Wald interval, double std dev (wider)"
df_params[!,:prior] .= "Wald interval, double std dev (wider)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

opts["pP0"] = Uniform(0, 0.1)
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Uniform(0, 0.1) (wider)"
df_params[!,:prior] .= "Uniform(0, 0.1) (wider)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

CSV.write("paper/outputs/s1-sensitivity/s1_sensitivity_pP0_states.csv", df_states_all)
CSV.write("paper/outputs/s1-sensitivity/s1_sensitivity_pP0_params.csv", df_params_all)


# ------------ TEST PRIOR DISTRIBUTION OF sigma ------------

# Fetch opts
opts = fetchDefaultOpts("binomial")

# Specify the initial distribution for Pt
phat = max(Y.nPos[1], 1)/Y.nSamples[1]
sd = sqrt((phat*(1-phat))/Y.nSamples[1])
opts["pP0"] = Truncated(Normal(phat, sd), 0, 1)

df_states_all = DataFrame()
df_params_all = DataFrame()

opts["paramPriors"] = [Uniform(0, 0.2)]
opts["initialParamSamplers"] = [Uniform(0, 0.2)]
opts["paramLimits"] = [(0, 0.2)]
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Uniform(0, 0.2) (default)"
df_params[!,:prior] .= "Uniform(0, 0.2) (default)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

opts["paramPriors"] = [Uniform(0, 1.0)]
opts["initialParamSamplers"] = [Uniform(0, 0.2)]
opts["paramLimits"] = [(0, 1.0)]
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "Uniform(0, 1.0) (wider)"
df_params[!,:prior] .= "Uniform(0, 1.0) (wider)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

opts["paramPriors"] = [Truncated(Cauchy(0, 100), 0, Inf)]
opts["initialParamSamplers"] = [Uniform(0.005, 0.02)]
opts["paramLimits"] = [(0, 100)]
(df_states, df_params, θ, diag) = fitModel(BinomialModel, Y, opts; skipResamplingPMMH=true, verbose=true, posteriorPredictive="binomial", checkLogLikStdDev=false)
df_states[!,:prior] .= "HalfCauchy(0, 100) (uninformative)"
df_params[!,:prior] .= "HalfCauchy(0, 100) (uninformative)"
df_states_all = vcat(df_states_all, df_states)
df_params_all = vcat(df_params_all, df_params)

CSV.write("paper/outputs/s1-sensitivity/s1_sensitivity_sigma_states.csv", df_states_all)
CSV.write("paper/outputs/s1-sensitivity/s1_sensitivity_sigma_params.csv", df_params_all)