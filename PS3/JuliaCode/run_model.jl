using Distributed, SharedArrays

#add processes
workers()
addprocs(2)


@Distributed.everywhere include("./PS3/JuliaCode/conesa_kueger.jl");

#prim, res = Initialize(); #=
#@time V_ret(prim, res);
#@time V_workers(prim, res);
#@time V_Fortran(res.r, res.w, res.b);
#@time SteadyStateDist(prim, res); =#

#agridf, consumption = V_Fortran(res.r, res.w, res.b);
#=
agridf
prim.a_grid
hcat(prim.a_grid, agridf, prim.a_grid - agridf)
=#
@time out_prim, out_res = MarketClearing(use_Fortran=false);

using Plots
theme(:juno)
plot(out_res.val_fun[:,:, end])
plot!(out_res.val_fun[:,:, end-1])
plot!(out_res.val_fun[:,:, end-2])



plot(out_prim.a_grid, out_res.val_fun[:, 1, end])
plot!(out_prim.a_grid, out_res.val_fun[:, 2, end])
plot(out_prim.a_grid, out_res.val_fun[:, 1, 20])
plot!(out_prim.a_grid, out_res.val_fun[:, 2,20])

plot(out_prim.a_grid, out_res.pol_fun[:,1,20])
plot!(out_prim.a_grid,out_res.pol_fun[:,2,20])



plot(out_prim.a_grid, out_res.val_fun[:, 1, 34])
plot!(out_prim.a_grid, out_res.val_fun[:, 2, 34])

# Calculate savings
savings = out_res.pol_fun[:,:,:] .- out_prim.a_grid;

plot(out_prim.a_grid, savings[:,1,20])
plot!(out_prim.a_grid, savings[:,2,20])

plot(out_prim.a_grid, out_res.pol_fun[:,1,20])
plot!(out_prim.a_grid, out_res.pol_fun[:,2,20])
plot!(out_prim.a_grid, out_prim.a_grid)

plot(savings[500,1,:])
plot!(savings[500,2,:])
plot!

# graph distributions

# accross assets for workers and retirees
a_dist = sum(out_res.F, dims = 2:3)
plot(out_prim.a_grid, a_dist[:, 1, 1])
