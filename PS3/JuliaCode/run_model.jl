using Plots, Distributed, SharedArrays

#add processes
#workers()
#addprocs(3)

theme(:juno)
@Distributed.everywhere include("PS3/JuliaCode/conesa_kueger.jl");

prim, res = Initialize();
#=
@elapsed V_ret(prim, res);
@elapsed V_workers(prim, res);
@elapsed SteadyStateDist(prim, res);
=#
@elapsed MarketClearing(prim, res);

plot(res.pol_fun[500,1,:])

plot(prim.a_grid, res.val_fun[:, 1, end])
plot!(prim.a_grid, res.val_fun[:, 2, end])
plot(prim.a_grid, res.val_fun[:, 1, 20])
plot!(prim.a_grid, res.val_fun[:, 2,20])

res.val_fun[:,:,end]

plot(prim.a_grid, res.val_fun[:, 1, 34])
plot!(prim.a_grid, res.val_fun[:, 2, 34])

# Calculate savings
savings = res.pol_fun[:,:,:] .- prim.a_grid;

plot(prim.a_grid, savings[:,1,20])
plot!(prim.a_grid, savings[:,2,20])

plot(prim.a_grid, res.pol_fun[:,1,20])
plot!(prim.a_grid, res.pol_fun[:,2,20])
plot!(prim.a_grid, prim.a_grid)

plot(savings[500,1,:])
plot!(savings[500,2,:])
plot!