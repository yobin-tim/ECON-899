using Plots
theme(:juno)
include("conesa_kueger.jl");

prim, res = Initialize();

V(prim, res);

plot(res.pol_fun[500,1,:])

plot(prim.a_grid, res.val_fun[:,1,50])
plot!(prim.a_grid, res.val_fun[:,2,50])

# Calculate savings
savings = res.pol_fun[:,:,:] .- prim.a_grid;

plot(prim.a_grid, savings[:,1,20])
plot!(prim.a_grid, savings[:,2,20])

plot(prim.a_grid, res.pol_fun[:,1,1])
plot!(prim.a_grid, res.pol_fun[:,2,1])