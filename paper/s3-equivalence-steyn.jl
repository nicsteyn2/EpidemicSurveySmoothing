
using CSV, DataFrames, Plots

include("../src/Simulations.jl")
include("../src/Opts.jl")
include("../src/Models.jl")
include("../src/FitModel.jl")
include("../src/LoadData.jl")

# Test logit model
(Y, tMax) = loadReactData(14, 19)

# Add steps at start and end
Y = vcat(DataFrame(t=[-1,-2,0], nSamples=[0,0,0], nPos=[0,0,0], phat=0, ssw=NaN, date=[Date("2021-09-06"), Date("2021-09-07"), Date("2021-09-08")], Round=NaN), Y)
Y = vcat(Y, DataFrame(t=[205,206,207], nSamples=[0,0,0], nPos=[0,0,0], phat=0, ssw=NaN, date=[Date("2022-04-01"), Date("2022-04-02"),Date("2022-04-03")], Round=NaN))

opts = fetchDefaultOpts("betabinom")
opts["chunkSize"] = 200
opts["pP0"] = Uniform(0.001, 0.2)
opts["paramPriors"] = [InverseGamma(0.0001, 0.0001), Uniform(0, 0.01)]
opts["initialParamSamplers"] = [Uniform(0.01, 0.018), Uniform(1e-4, 3e-4)]

st = time()
(df_states, df_params, θ, diag) = fitModel(logitModel, Y, opts; skipResamplingPMMH=true, verbose=true, checkLogLikStdDev=false)
en = time()
println("Time taken: ", en-st)

Chains(θ*1e4)
