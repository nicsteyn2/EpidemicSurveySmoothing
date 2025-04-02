
using Dates, CSV, DataFrames


function loadReactData(min_round, max_round)

    # Load the file, correctly set the dates, filter by round, and set number of days since first obs
    df_raw = DataFrame(CSV.File("data/reactdata.csv"))
    df_raw = df_raw[(df_raw.round.>=min_round) .& (df_raw.round.<=max_round),:]
    df_raw.t = Dates.value.(df_raw.date - minimum(df_raw.date)) .+ 1

    # Get the maximum time
    tMax = maximum(df_raw.t)

    # Now create the output file (we have to fill in the zeros)
    t = collect(1:tMax)
    dates = minimum(df_raw.date) .+ Day.(t .- 1)

    nSamplesNew = Int.(zeros(tMax))
    nPositiveNew = Int.(zeros(tMax))
    phatNew = zeros(tMax)
    sswNew = zeros(tMax)
    Round = zeros(tMax) .+ NaN

    nSamplesNew[df_raw.t] = Int.(df_raw.nSamples)
    nPositiveNew[df_raw.t] = Int.(df_raw.nPos)
    phatNew[df_raw.t] = df_raw.phat
    sswNew[df_raw.t] = df_raw.ssw
    Round[df_raw.t] = df_raw.round

    df_all = DataFrame(t=t, nSamples=nSamplesNew, nPos=nPositiveNew, phat=phatNew, ssw=sswNew, date=dates, Round=Round)

    return(df_all, tMax)

end


function loadVariantData(min_round, max_round)

    # Load the file, correctly set the dates, filter by round, and set number of days since first obs
    df_raw = DataFrame(CSV.File("data/reactdata_variants.csv"))
    df_raw = df_raw[(df_raw.round.>=min_round) .& (df_raw.round.<=max_round),:]
    df_raw.t = Dates.value.(df_raw.date - minimum(df_raw.date)) .+ 1

    # Get the maximum time
    tMax = maximum(df_raw.t)

    # Now create the output file (we have to fill in the zeros)
    t = collect(1:tMax)
    dates = minimum(df_raw.date) .+ Day.(t .- 1)

    nSamplesNew = zeros(tMax)
    nPositiveNew = zeros(tMax)
    nSequenceNew = zeros(tMax)
    nWildtypeNew = zeros(tMax)
    nAlphaNew = zeros(tMax)
    nDeltaNew = zeros(tMax)
    nOmicronNew = zeros(tMax)
    nOmicronBA2New = zeros(tMax)
    Round = zeros(tMax) .+ NaN

    nSamplesNew[df_raw.t] = df_raw.n
    nPositiveNew[df_raw.t] = df_raw.nPos
    nSequenceNew[df_raw.t] = df_raw.nSeq .- df_raw.nOther
    nWildtypeNew[df_raw.t] = df_raw.nWT
    nAlphaNew[df_raw.t] = df_raw.nAL
    nDeltaNew[df_raw.t] = df_raw.nDE
    nOmicronNew[df_raw.t] = df_raw.nOM
    nOmicronBA2New[df_raw.t] = df_raw.nBA2 .+ df_raw.nBA3
    Round[df_raw.t] = df_raw.round

    nSamplesNew = Int.(nSamplesNew)
    nPositiveNew = Int.(nPositiveNew)
    nSequenceNew = Int.(nSequenceNew)
    nWildtypeNew = Int.(nWildtypeNew)
    nAlphaNew = Int.(nAlphaNew)
    nDeltaNew = Int.(nDeltaNew)
    nOmicronNew = Int.(nOmicronNew)
    nOmicronBA2New = Int.(nOmicronBA2New)

    df_all = DataFrame(t=t, nSamples=nSamplesNew, nPos=nPositiveNew, nSeq=nSequenceNew, nWT=nWildtypeNew, nAL=nAlphaNew, nDE=nDeltaNew, nOM=nOmicronNew, nBA=nOmicronBA2New, date=dates, Round=Round)

    return(df_all, tMax)

end
