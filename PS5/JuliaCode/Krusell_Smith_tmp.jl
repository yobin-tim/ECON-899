@everywhere using Distributed, SharedArrays

workers()
addprocs(1)

@Distributed.everywhere include("./HelpfulFunctions.jl")

Epsilon_idio, Epsilon_agg = draw_shocks(shock, prim.N, prim.T)

