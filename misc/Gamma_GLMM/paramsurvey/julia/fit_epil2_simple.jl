using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "fitlib.jl"))

form = @formula(y ~ 1 + trt + (1 | subject))
run_survey("epil2_simple", form, Gamma(), LogLink(), :subject, 2)
