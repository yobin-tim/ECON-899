include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

market_clearing(prim, res)

using Plots
theme(:juno)

plot(res.W_val)
scatter!(res.W_val)
plot(res.x_opt)
scatter!(res.x_opt)

res.μ = Tμ(prim, res ) 
        
# Calculate optimal labor demand for each firm (for each productivity level)
n_opt = prim.n_optim.(prim.s_vals, res.p)
# Calculate profit for each firm (for each productivity level)
prof = prim.Π.(prim.s_vals, res.p, n_opt)
# Calculate  mass of firms in the market (for each productivity level)
mass = res.μ + res.M * prim.ν

# Calculate Total labor demand
tot_labor_demand = n_opt' * mass
# Calculate total profits 
tot_profit = prof' * mass
# Calculate total supply of labor
tot_labor_supply = 1/prim.A - tot_profit

err = tot_labor_demand - tot_labor_supply #reset error level

res.M =2.5

res.M = 0.75 * res.M


res.p = 0.7