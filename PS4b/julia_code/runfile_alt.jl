include("functions_alt.jl")

out_res=SolveModel()

#Problems 1 and 2
using Latexify, Plots
    println("Table for Problem 1:\n")
    Table=convert(DataFrame,
        hcat(0:35,out_res.EV_P1,out_res.EV_CCP))
        rename!(Table, [:state_id,:EV,:EV_CPP])
        Table[!,:state_id] = convert.(Int64,Table[!,:state_id])
    latexify(Table, env=:tabular, fmt = x -> round(x, digits=3)) |> print

#Problem 4
println("The optimal λ is $(out_res.λOpt)")
