using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "fitlib.jl"))

form = @formula(imps79 ~ 1 + TxDrug * Week + (1 | id))
run_survey("schizophrenia", form, Gamma(), LogLink(), :id, 4)
