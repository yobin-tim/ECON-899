using Parameters, LinearAlgebra

# Structure that holds the parameters for the model
@with_kw struct Primitives
    β           ::Float64          = 0.8
    θ           ::Float64          = 0.64
    s_vals      ::Array{Float64}   = [3.98e-4, 3.58, 6.82, 12.18, 18.79]
    nS          ::Int64            = length(s_vals)
    x_vals      ::Array{Int64}     = [0, 1]
    emp_lev     ::Array{Float64}   = [1.3e-9, 10, 60, 300, 1000]
    trans_mat   ::Array{Float64,2} =[0.6598  0.2600  0.0416  0.0331  0.0055 ;
                                    0.1997  0.7201  0.0420  0.0326  0.0056 ;
                                    0.2000  0.2000  0.5555  0.0344  0.0101 ;
                                    0.2000  0.2000  0.2502  0.3397  0.0101 ;
                                    0.2000  0.2000  0.2500  0.3400  0.0100]
    ν           ::Array{Float64}   = [0.37, 0.4631, 0.1102, 0.0504, 0.0063]
    A           ::Float64          = 1/200
    c_f         ::Int64            = 10
    c_e         ::Int64            = 5

    # Price grid
    p_min       ::Float64 = 0.01
    p_max       ::Float64 = 10
    # nP          ::Int64   = 10
    # p_grid      ::Array{Float64}   = range(p_min, stop = p_max, length = nP)

    # Optimal decision rules
    n_optim     ::Function         = ( s , p ) -> (θ * p * s) ^ (1/(1 - θ))

    Π           ::Function         = ( s , p,  n ) -> p*s*( n )^θ - n  - p * c_f

    # Limits for the mass of entrants
    M_min       ::Float64 = 1.0
    M_max       ::Float64 = 10.0


end

# Structure that stores results
mutable struct Results
    W_val   ::Array{Float64}         # Firm value given state variables
    n_opt   ::Array{Float64}         # Optimal labor demand for each possible state
    x_opt   ::Array{Int64}           # Optimal firm decicion for each state
    p       ::Float64                # Market clearing price
    μ       ::Array{Float64}         # Distribution of Firms
    M       ::Float64                # Mass of entrants
end

# Initialize model
function Initialize()
    prim  = Primitives()

    W_val = zeros(prim.nS)
    n_opt = zeros(prim.nS)
    x_opt = zeros(prim.nS)
    p = 7.0
    μ = ones(prim.nS) / prim.nS # Uniform distribution is the initial guess
    M = 5.0

    res = Results(W_val, n_opt, x_opt, p, μ, M)

    return prim, res

end

# Bellman operator for W
function W(prim::Primitives, res::Results)
    @unpack Π, n_optim,s_vals, nS, trans_mat, c_f, β = prim
    @unpack p = res

    temp_val = zeros(size(res.W_val))

    n_opt = prim.n_optim.(s_vals, p)
    profit_state = Π.(s_vals, p, n_opt)

    # Iterate over all possible states
    for s_i ∈ 1:nS

        prof = profit_state[s_i]

        # Calculate expected continuation value

        exp_cont_value = trans_mat[s_i, :]' * res.W_val

        x = ( exp_cont_value > 0 ) ? 0 : 1

        temp_val[s_i] = prof + β * (1 - x) * (exp_cont_value )
        res.x_opt[s_i] = x

    end

    res.W_val = temp_val

end # W

#Value function iteration for W operator
function TW_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4)

    n = 0 #counter
    err = 100.0 #initialize error
    while  (err > tol) & (n < 4000)#begin iteration
        W_val_old = copy(res.W_val)
        W(prim, res)
        err = maximum(  abs.(W_val_old - res.W_val ) ) #reset error level
        n+=1
        if n % 100 == 0
            println("Iter =", n , " Error = ", err)
        end
    end
end # TW_iterate

function market_clearing(prim::Primitives, res::Results; tol::Float64 = 1e-3,  n_max::Int64 = 1000)
    @unpack Π, nS, trans_mat, ν, c_e, p_min, p_max = prim

    θ = 0.99
    n = 0

    while n < n_max
        TW_iterate(prim, res)
        # W(prim, res)
        # Calculate EC
        EC = sum(res.W_val .* ν)/res.p - c_e  #This was previously - res.p * c_e, but that doesn't seem right?
        # println("p = ", res.p," EC = ", EC, " tol = ", tol)
        if abs(EC) > tol * 10000
        # adjust tuning parameter based on EC
            θ = 0.5
        elseif abs(EC) > tol * 5000
            θ = 0.75
        elseif abs(EC) > tol * 1000
            θ = 0.9
        else
            θ = 0.99
        end
        # println(n+1, " iterations; EC = ", EC, ", p = ", res.p, ", p_min = ", p_min, ", p_max = ", p_max, ", θ = ", θ)
        if abs( EC ) < tol
            println("Market Cleared in $(n+1) iterations.")
            break
        end

        # adjust price toward bounds according to tuning parameter
        if EC > 0
            p_old = res.p
            res.p = θ*res.p + (1-θ)*p_min
            p_max = p_old
        else
            p_old = res.p
            res.p = θ*res.p + (1-θ)*p_max
            p_min = p_old
        end

        n += 1

    end
end

# Calculate distribution of firms
function Tμ(prim::Primitives, res::Results)
    @unpack trans_mat, nS, ν = prim

    stay = 1 .- res.x_opt
    μ₁ = zeros(nS)
    for s_i ∈ 1:nS
        # Add the mass of firms that stay in the market
        μ₁[s_i] = trans_mat[s_i, :]' * (stay .* res.μ)
        # and the mass of firms that enter the market
        μ₁[s_i] += trans_mat[s_i, :]' * (stay .*  ν * res.M)
    end
    return μ₁
end #

# Calculate steady state distribution of firms
# TODO: Check why this is not necesary
function Tμ_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-3, n_max::Int64 = 1000)

    n = 0 #counter
    err = 100.0 #initialize error
    while  (err > tol) & (n < n_max)#begin iteration
        μ_new = Tμ(prim, res)
        err = maximum( abs.( μ_new .- res.μ ) ) #reset error level
        res.μ = μ_new
        n+=1
        if n % 100 == 0
            println("Iter =", n , " Error = ", err)
        end
        if n == 4000
            println("Iteration limit reached")
        end
        if err < tol
            println("Distribution converged in $n iterations")
            break
        end
    end

end # Tμ_iterate

# Iterate until labor market clears
function Tμ_iterate_until_cleared(prim::Primitives, res::Results;  tol::Float64 = 1e-3, n_max::Int64 = 1000)
    @unpack Π, n_optim , s_vals, ν, A, M_min, M_max = prim

    θ = 0.5
    n = 0 #counter
    err = 100.0 #initialize error

    # Find Z Matrix
    Z = repeat((1 .- res.x_opt)', prim.nS)

    # Calculate optimal labor demand for each firm (for each productivity level)
    n_opt = prim.n_optim.(prim.s_vals, res.p)
    # Calculate profit for each firm (for each productivity level)
    prof =  prim.Π.( prim.s_vals, res.p, n_opt)

    while  (abs(err) > tol) & (n < n_max)#begin iteration

        res.μ = res.M *((Z - I)^(-1)) * Z*ν # Find the distribution of firms using initial guess for entrants

        # Calculate  mass of firms in the market (for each productivity level)
        mass = res.μ + res.M * ν

        # Calculate Total labor demand
        tot_labor_demand = n_opt' * mass
        # Calculate total profits
        tot_profit = prof' * mass
        # Calculate total supply of labor
        tot_labor_supply = 1/A - tot_profit

        # Labor MarketClearing condition
        LMC = tot_labor_demand - tot_labor_supply

        if abs(LMC) > tol * 10000
            # adjust tuning parameter based on LMC
                θ = 0.5
            elseif abs(LMC) > tol * 5000
                θ = 0.75
            elseif abs(LMC) > tol * 1000
                θ = 0.9
            else
                θ = 0.99
            end

        println(n+1, " iterations; LMC = ", LMC, ", M = ", res.M, ", M_min = ", M_min, ", M_max = ", M_max, ", θ = ", θ)

        if abs( LMC ) < tol
            println("Labor Market Cleared")
            break
        end
        # adjust price toward bounds according to tuning parameter
        if LMC > 0
            M_old = res.M
            res.M = θ*res.M + (1-θ)*M_min
            M_max = M_old
        else
            M_old = res.M
            res.M = θ*res.M + (1-θ)*M_max
            M_min = M_old
        end

        n += 1
    end

end
