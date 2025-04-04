
include("../src/Simulations.jl")
include("../src/Opts.jl")
include("../src/Models.jl")
include("../src/FitModel.jl")
include("../src/LoadData.jl")

# # Variants 8 to 13
# (Y, tMax) = loadVariantData(8, 13)
# θ = [0.03, 0.01, 0.012, 1.5e-4]
# opts = fetchDefaultOpts("variants8to13")
# opts["N"] = 4000
# opts["posteriorNumberOfParticles"] = 10000
# opts["chunkSize"] = 200

# function fetchInitialPrevalenceDist(nPos, nSamples)
#     phat = max(nPos, 1)/nSamples
#     sd = sqrt((phat*(1-phat))/nSamples)
#     return(Truncated(Normal(phat, sd), 0, 1.0))
# end

# opts["pP0-WT"] = fetchInitialPrevalenceDist(Y.nWT[4], Y.nSamples[4])
# opts["pP0-AL"] = fetchInitialPrevalenceDist(Y.nAL[4], Y.nSamples[4])
# opts["pP0-DE"] = fetchInitialPrevalenceDist(Y.nDE[119], Y.nSamples[119])

# # test = VariantsModel(θ, Y, opts)
# # (sd, ll) = estimateStdDevLogLik(100, VariantsModel8to13, θ, Y, opts; showProgress=true, ignoreerror=false)
# (df_states, df_params, θ, diag) = fitModel(VariantsModel8to13, Y, opts; skipResamplingPMMH=true, verbose=true, checkLogLikStdDev=false)

# CSV.write("paper/outputs/s10-variants/variants_r8to13_results.csv", df_states)
# CSV.write("paper/outputs/s10-variants/variants_r8to13_params.csv", df_params)


# Variants 14 to 19
(Y, tMax) = loadVariantData(14, 19)
θ = [0.015, 0.035, 0.015, 2e-4]
opts = fetchDefaultOpts("variants14to19")
opts["N"] = 10000
opts["posteriorNumberOfParticles"] = 10000
opts["chunkSize"] = 200

opts["pr0-DE"] = Uniform(-0.1, 0.1)
opts["pr0-OM"] = Uniform(0, 0.4)
opts["pr0-BA"] = Uniform(0, 0.4)

opts["pP0-DE"] = fetchInitialPrevalenceDist(Y.nDE[1], Y.nSamples[1])
opts["pP0-OM"] = fetchInitialPrevalenceDist(1, 100)
opts["pP0-BA"] = fetchInitialPrevalenceDist(1, 100)

# test = VariantsModel14to19(θ, Y, opts)
# (sd, ll) = estimateStdDevLogLik(100, VariantsModel14to19, θ, Y, opts; showProgress=true, ignoreerror=true)
(df_states, df_params, θ, diag) = fitModel(VariantsModel14to19, Y, opts; skipResamplingPMMH=true, verbose=true, checkLogLikStdDev=false)

CSV.write("paper/outputs/s10-variants/variants_r14to19_results.csv", df_states)
CSV.write("paper/outputs/s10-variants/variants_r14to19_params.csv", df_params)