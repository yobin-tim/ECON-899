prim, res, shocks = Initialize()

# Here I'm running an iteration by hand to see what happens
zi = 1
z = prim.z_val[zi]   # productivity
L = (1 - prim.u[zi])*prim.ē  # aggregate effective labor
Ki = 1
K = prim.K_grid[Ki]
r, w = prim.r_mkt(K, L, z), prim.w_mkt(K, L, z)
Knext = prim.k_forecast(z, res.a, res.b, K)
ei = 1
cand_last = 1
ezi = 2*(zi - 1) + ei
e = prim.e_val[ei]       # employment status
ki = 1
k = prim.k_grid[ki]
cand_max = -Inf     # 
pol_max = 0         # policy function maximum
budget   = r*k + w*e*prim.ē + (1-prim.δ)*k
#run the fisrt iteration of the innermost loop
kpi = cand_last
knext = prim.k_grid[kpi]
c = budget - knext
# The important part is to use the interpolated functions

# Try the Bellman equation
@time Bellman(prim, res, shocks)

# the same for V_iterate
err = 100
tol = 1e-3
n = 1
v_next = copy(res.val_fun)
Bellman(prim, res, shocks)
maximum(abs.(v_next .- res.val_fun))

V_iterate(prim, res, shocks)

using Plots
theme(:juno)
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style
plot(prim.k_grid, res.val_fun[:,1,1,1])
plot!(prim.k_grid, res.val_fun[:,1,1,2])
plot(prim.k_grid, res.val_fun[:,1,2,1])
plot!(prim.k_grid, res.val_fun[:,1,2,2])
