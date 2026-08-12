using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "fitlib.jl"))

form = @formula(crate ~ 1 + (1 | location) + (1 | fyear))
run_survey("report4bb", form, Gamma(), LogLink(), (:location, :fyear), 1)
