include("./krusell_smith.jl")
prim, res, shocks = Initialize()

Bellman(prim, res, shocks)
p1 = copy(res.pol_fun)

res.pol_fun[:,:,1,1]
p1[:,:,1,1]

res.k_forecast_grid = 10.0 * ones(prim.nK, prim.nZ)

V_iterate(prim, res, shocks; use_dense_grid=true)

V1 = Simulation(prim, res, shocks)

reg_coefs, R_squared = auto_reg(prim, res, shocks)

err = sum( (res.a - reg_coefs[1]).^2) + sum( (res.b - reg_coefs[2]).^2 )

# Update regression coefficients
λ = 0.5
res.a = reg_coefs[1] * (1 - λ) + λ * res.a  
res.b = reg_coefs[2] * (1 - λ) + λ * res.b

k_forecast_grid_old = copy(res.k_forecast_grid)

# Re-calculate forcasted capital values
k_forecast_grid = zeros(prim.nK, prim.nZ)
k_forecast_grid[:, 1] = prim.k_forecast.(1, Ref(res.a), Ref(res.b), prim.K_grid)
k_forecast_grid[:, 2] = prim.k_forecast.(2, Ref(res.a), Ref(res.b), prim.K_grid)
res.k_forecast_grid = copy(k_forecast_grid)

pol_old = copy(res.pol_fun)

res.k_forecast_grid = 14 * ones(prim.nK, prim.nZ)

V_iterate(prim, res, shocks; use_dense_grid=true)

p2 = res.pol_fun[:,:,1,1]

res.val_fun_interp[(1,1)](10.45, 11)
res.val_fun_interp[(1,1)](10.45, 14)

res.val_fun_interp[(1,1)](10.45, prim.K_min:0.01:prim.K_max)



V2 = Simulation(prim, res, shocks)



auto_reg(prim, res, shocks)

using Plots
theme(:juno)
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

 plot(res.pol_fun[:,10,1,1])
plot!(res.pol_fun[:,10,1,2])
plot!(res.pol_fun[:,10,2,1])
plot!(res.pol_fun[:,10,2,2])


k_grid_dense = range(prim.k_min, stop=prim.k_max, length=2000)

plot(res.pol_fun_interp[1,1].(k_grid_dense, prim.K_grid[1]))


SolveModel(1)

Bellman(prim, res, shocks)

@unpack k_grid, K_grid, nk, nK, nZ, nE, ē, w_mkt, r_mkt, β, δ, k_forecast, z_val, e_val, u, y, util, k_min, k_max = prim
@unpack a, b, val_fun, val_fun_interp, k_forecast_grid = res
@unpack Π = shocks

zi = 1
# save aggregate shock and relevant variables
z = z_val[zi]   # productivity
L = (1 - u[zi])*ē         # aggregate effective labor
Ki = nK
K = K_grid[Ki]  # aggregate capital
r, w = r_mkt(K, L, z), w_mkt(K, L, z)
# estimate next period capital 
# Knext = k_forecast(z, a, b, K)
Knext = k_forecast_grid[Ki, zi]
Knext = min(Knext, prim.K_max)
ei = 1
ezi = 2*(zi - 1) + ei
e = e_val[ei]       # employment status
ki = nk
k = k_grid[ki]      # current period capital
cand_max = -Inf     # intial value maximum
pol_max  = 1        # policy function maximum
budget   = r*k + w*e*ē + (1-δ)*k
k_grid_dense = range(prim.k_min, stop= prim.k_max, step=0.01)
c_pos = budget .- k_grid_dense
# k_grid_dense = k_grid_dense[c_pos .> 0 ]
# c_pos = c_pos[c_pos .> 0 ]
utils = util.( c_pos )
fut_vals = [res.val_fun_interp[(i, j)].(k_grid_dense, Knext) for i ∈ 1:2 for j ∈ 1:2];
fut_vals = hcat(fut_vals...)
exp_val_next = [shocks.Π[ezi, :]' * fut_vals[i,:] for i ∈ 1:size(fut_vals,1)] 
cont_val = utils + β*exp_val_next

plot(k_grid_dense, cont_val)

cand_last = k_grid_dense[argmax(cont_val)]
cand_max = maximum(cont_val)
pol_max  = argmax(cont_val)