include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

market_clearing(prim, res; n_max =100)

# Testing Alex Idea
# Calculate A Matrix 
using LinearAlgebra
y = 1 .- res.x_opt
A = repeat(y', prim.nS) .* prim.trans_mat'
Tμ(M) = M *((A - I)^(-1)) * A*prim.ν

μ = Tμ(1)
n_opt = prim.n_optim.(prim.s_vals, res.p)
prof =  prim.Π.( prim.s_vals, res.p, n_opt)
mass = μ + 1 *  prim.ν

tot_labor_demand = n_opt' * mass 

sum([prim.n_optim.(prim.s_vals[i], res.p)*(μ[i] + 1*prim.ν[i]) for i∈1:prim.nS])

# Calculate total profits 
tot_profit = prof' * mass
# Calculate total supply of labor
tot_labor_supply = 1/prim.A - tot_profit


function labor_supply_demand(prim::Primitives, res::Results, M::Float64)
    @unpack c_e, ν = prim

    μ = Tμ(M)
    
    # Calculate optimal labor demand for each firm (for each productivity level)
    n_opt = prim.n_optim.(prim.s_vals, res.p)
    # Calculate profit for each firm (for each productivity level)
    prof =  prim.Π.( prim.s_vals, res.p, n_opt)
    # Calculate  mass of firms in the market (for each productivity level)
    mass = μ + M *  prim.ν

    # Calculate Total labor demand
    tot_labor_demand = n_opt' * mass 

    # Calculate total profits 
    tot_profit = prof' * mass
    # Calculate total supply of labor
    tot_labor_supply = 1/prim.A - tot_profit

    return tot_labor_supply, tot_labor_demand

end

M_range = range(prim.M_min, stop = 5, length=100)

L_S_D = labor_supply_demand.(Ref(prim), Ref(res), M_range)
L_S_D = hcat(collect.(L_S_D)...)

plot(M_range, L_S_D[1,:], c=1, label = "Labor Supply", legend=:topleft)
plot!(M_range, L_S_D[2,:], c=2, label = "Labor Demand")

function calculate_EC(prim::Primitives, res::Results, p::Float64)
    @unpack c_e, ν = prim
    res.p = p
    TW_iterate(prim, res)
    # Calculate EC
    EC = sum(res.W_val .* ν) - res.p * c_e
    return EC
end

p_grid = prim.p_min:0.5:prim.p_max 
ECs = calculate_EC.(Ref(prim), Ref(res), p_grid)

using Interpolations, Optim

itp = CubicSplineInterpolation(M_range, abs.(L_S_D[1,:] - L_S_D[2,:]))

plot(M_range, itp.(M_range))

opt = optimize( x -> abs(itp(x)),  M_range[1], M_range[end] )



using Plots
theme(:juno)

plot(prim.s_vals, res.W_val)
scatter!(prim.s_vals, res.W_val, c=1)
plot(prim.s_vals, res.x_opt)
scatter!(prim.s_vals, res.x_opt, c=1)
xticks!(prim.s_vals)
res.M = 3