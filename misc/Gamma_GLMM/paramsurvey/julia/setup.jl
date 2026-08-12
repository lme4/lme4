using Pkg
Pkg.activate(@__DIR__)
Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaStats/MixedModels.jl", rev="pa/dispersion-again"))
Pkg.add(["CSV", "DataFrames", "GLM", "StatsModels", "DaemonMode"])
Pkg.instantiate()
Pkg.precompile()
using MixedModels, CSV, DataFrames, GLM, StatsModels, DaemonMode
println("setup OK, MixedModels loaded from: ", pathof(MixedModels))
