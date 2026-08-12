using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "fitlib.jl"))

form = @formula(y ~ 1 + Base * trt + Age + Visit + (1 + Visit | subject))
run_survey("epil2_complex", form, Gamma(), LogLink(), :subject, 6)
