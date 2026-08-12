## Shared fitting/extraction routine for the paramsurvey Julia arm
## (MixedModels.jl, pa/dispersion-again). Mirrors the R side's
## resultsToDF() column shape (i, status, singular, msg, time_sec, sd1,
## sd2, corr, phi, negll, beta1..betaP) so 10_ingest_julia_results.R can
## read any dataset's output the same way.

using MixedModels, CSV, DataFrames, GLM, StatsModels

## sdcorr_spec: either a single Symbol (one grouping factor -- 1 or 2
## correlated terms, extracted like R's extract_sdcorr()) or a
## Tuple{Symbol,Symbol} (two independent single-term grouping factors,
## no correlation -- matches Report4BB's (1|location)+(1|fyear) shape).
function extract_sdcorr(gm, sdcorr_spec)
    vc = VarCorr(gm).σρ
    if sdcorr_spec isa Tuple
        g1, g2 = sdcorr_spec
        return (sd1 = only(getproperty(vc, g1).σ), sd2 = only(getproperty(vc, g2).σ),
                corr = NaN)
    else
        g = getproperty(vc, sdcorr_spec)
        sd = g.σ
        if length(sd) == 1
            return (sd1 = sd[1], sd2 = NaN, corr = NaN)
        else
            return (sd1 = sd[1], sd2 = sd[2], corr = only(g.ρ))
        end
    end
end

function run_survey(example::String, form, family, link, sdcorr_spec, nbeta::Int;
                     datadir=joinpath(@__DIR__, "data", example),
                     outfile=joinpath(@__DIR__, "results_$(example).csv"))
    reps = sort(filter(f -> endswith(f, ".csv"), readdir(datadir)))
    B = length(reps)
    println("Fitting MixedModels.jl ($family, $link) to $B replicates of $example")

    rows = Vector{NamedTuple}(undef, B)
    for (i, fname) in enumerate(reps)
        dat = CSV.read(joinpath(datadir, fname), DataFrame)
        t0 = time()
        status = "error"; msg = ""; singular = false
        sd1 = sd2 = corr = phi = negll = NaN
        betas = fill(NaN, nbeta)
        try
            gm = fit(MixedModel, form, dat, family, link; progress=false)
            rv = gm.optsum.returnvalue
            status = rv in (:FTOL_REACHED, :SUCCESS, :XTOL_REACHED, :ROUNDOFF_LIMITED) ? "clean" : "warning"
            msg = string(rv)
            singular = issingular(gm)
            sc = extract_sdcorr(gm, sdcorr_spec)
            sd1, sd2, corr = sc.sd1, sc.sd2, sc.corr
            phi = dispersion(gm, true)
            negll = deviance(gm)
            betas = collect(gm.β)
        catch e
            status = "error"
            msg = sprint(showerror, e)
        end
        time_sec = time() - t0
        rows[i] = merge((i = i, status = status, singular = singular, msg = msg,
                          time_sec = time_sec, sd1 = sd1, sd2 = sd2, corr = corr,
                          phi = phi, negll = negll),
                         NamedTuple{Tuple(Symbol.("beta", 1:nbeta))}(Tuple(betas)))
        println("  rep $i: status=$status singular=$singular time=$(round(time_sec, digits=3))s")
    end

    df = DataFrame(rows)
    CSV.write(outfile, df)
    println("saved to $outfile")
    println(combine(groupby(df, :status), nrow))
    df
end
