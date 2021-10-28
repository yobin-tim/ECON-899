include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

@time begin
market_clearing(prim, res; n_max =100)
println(res.p)
end


# Testing an alternative approach to the market clearing problem
# I will evaluate in a few point and then interpolate to find the optimal price
using Interpolations, Optim

function calculate_EC(prim::Primitives, res::Results, price::Float64)
    @unpack c_e, ν = prim
    res.p = price
    TW_iterate(prim, res)
    EC = sum(res.W_val .* ν)/ res.p -c_e
    return EC
end

function find_price(prim::Primitives, res::Results, n_points::Int64)
    @unpack c_e, ν, p_min, p_max = prim
    # Define the range of prices to be evaluated
    price_range = range(p_min, stop = p_max, length=n_points)
    # Calculate the EC for each price
    ECs = calculate_EC.(Ref(prim), Ref(res), price_range)
    # Interpolate to find the objective function
    itp = CubicSplineInterpolation(price_range, ECs)
    # Find the minimum of the objective function
    opt = optimize( x -> abs(itp(x)),  price_range[1], price_range[end] )
    return opt.minimizer
end

@time find_price(prim, res, 1000)

using BenchmarkTools

@benchmark begin
    prim, res = Initialize()
    market_clearing(prim, res; n_max =100)
end

@benchmark begin
    prim, res = Initialize()
    find_price(prim, res, 1000)
end

# Conclusion: For this problem, the alternative approach is much slower.

prim, res = Initialize()
market_clearing(prim, res; n_max =100)
Tμ_iterate_until_cleared(prim, res; n_max =500)

A = (1 .- res.x) 

function labor_supply_demand(prim::Primitives, res::Results, M::Float64)
    @unpack c_e, ν = prim
    res.M = M
    

    
    # Calculate optimal labor demand for each firm (for each productivity level)
    n_opt = prim.n_optim.(prim.s_vals, res.p)
    # Calculate profit for each firm (for each productivity level)
    prof =  prim.Π.( prim.s_vals, res.p, n_opt)
    # Calculate  mass of firms in the market (for each productivity level)
    mass = res.μ + res.M *  prim.ν
    mass = mass/sum(mass)

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

ECs = calculate_EC.(Ref(prim), Ref(res), price_range)

itp = CubicSplineInterpolation(price_range, ECs)

plot(price_range, itp.(price_range) .+ prim.c_e)
plot!([price_range[1], price_range[end]], [prim.c_e, prim.c_e], c=2)
xticks!(0.2:0.1:1.3)

plot(price_range, abs.(itp.(price_range)))

opt = optimize( x -> abs(itp(x)),  price_range[1], price_range[end] )



using Plots
theme(:juno)

plot(prim.s_vals, res.W_val)
scatter!(prim.s_vals, res.W_val, c=1)
plot(prim.s_vals, res.x_opt)
scatter!(prim.s_vals, res.x_opt, c=1)
xticks!(prim.s_vals)
res.M = 3