using Distributed, SharedArrays, JLD

#add processes
workers()
addprocs(2)

@Distributed.everywhere include("./conesa_kueger_tmp.jl");

prim, res = Initialize()

out_prim, out_res = ShootBackward(prim, res)
## Exercise 1 1.
##out_prim, out_res = MarketClearing(use_Fortran=false, tol = 1e-3);
# prim_noSS, res_noSS = MarketClearing(use_Fortran=false, tol = 1e-3, ss = false);

# save("../Data/Initial_Conditions.jld",
#      "Γ_0", out_res.F,
#      "V_0", out_res.val_fun,
#      "K_θ", out_res.K,
#      "K", res_noSS.K,
#      "L_θ", out_res.L,
#      "L", res_noSS.L)



