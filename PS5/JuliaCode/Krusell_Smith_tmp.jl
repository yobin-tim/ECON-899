using Distributed, SharedArrays, Interpolations, Optim

workers()
addprocs(1)

@Distributed.everywhere include("./HelpfulFunctions.jl")

prim, grid, shock, res = Initialize();

ℇ_idio, ℇ_agg = draw_shocks(shock, prim.N, prim.T)

pf_k_up, pf_v_up = Bellman(prim, grid, shock, res)
