<<<<<<< HEAD
#########################################
# Problem Set 2, ECON 899
# Goal Due Date: 9/20/2021
#########################################

#ReadMe: This file contains functions and variables. Run the code from the file PS2 compute.
# using Plots, Parameters, Pkg
#Keyword-enabled structure to hold model primitives

@with_kw struct Primitives
    β::Float64 = 0.9932 #discount factor
    α::Float64 = 1.5 #RRA coefficient
    e::Float64 = 1 #earnings, when employed
    u::Float64 = 0.5 #earnings, when unemployed
    nS::Int64 = 2 #number of states
    Πee::Float64 = 0.97
    Πuu::Float64 = 0.5
    Πue::Float64 = 0.03
    Πeu::Float64 = 0.5
    Π::Array{Float64, 2} = collect([Πee Πeu; Πue Πuu]) #in case I choose matrix indexing
    nA::Int64 = 1000 #number of asset grid points
    A_min::Int64 = -2 #asset lower bound
    A_max::Int64 = 5 #asset upper bound
    A_grid::Array{Float64, 1} = collect(range(A_min, length = nA, stop = A_max)) #asset grid
    a_lbar::Float64 = -2
    a_ubar::Float64 = 30
end

#Structure that holds model results
mutable struct Results
    val_func::Array{Float64, 2} #value function
    pol_func::Array{Float64, 2} #policy function
    q::Float64
    μ::Array{Float64,2}
end


#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.nS, prim.nA) #initial value function guess
    pol_func = zeros(prim.nS, prim.nA) #initial policy function guess
    q = (1 + prim.β)/2 #initial bondprice guess
    μ = zeros(prim.nS, prim.nA)
    for i = 1:prim.nS
        for j = 1:prim.nA
            μ[i,j] = (prim.A_grid[j] - prim.a_lbar)/(prim.a_ubar - prim.a_lbar)*prim.Π[i,i]
        end
    end

    res = Results(val_func, pol_func, q, μ) #initialize results struct
    prim, res #return deliverables
end

## Function for solving for the decision rule
function Bellman(prim::Primitives,res::Results)
    @unpack val_func, q = res
    @unpack β, α, e, u, Π, A_grid, nS, nA = prim
    v_next = zeros(nS, nA)

    for S_index = 1:nS
        if S_index == 1
            S = e
        else S = u
        end
        for A_index = 1:nA
            A = A_grid[A_index]
            candidate_max = -Inf
            budget = S + A
            for i = 1:nA
                c = budget - q * A_grid[i]
                if c > 0
                    val = (c^(1-α)-1)/(1 - α)
                else
                    val = -Inf
                end
                    for j = 1:nS
                        val = val + β * Π[j, S_index] * val_func[j, i]
                    end
                    if val > candidate_max
                        candidate_max = val
                        res.pol_func[S_index, A_index] = A_grid[i]
                    end
            end
            v_next[S_index, A_index] = candidate_max
        end
    end
    v_next
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter

    while err>tol && n <= 1000 #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func)) #reset error level
        res.val_func = v_next #update value function
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

## Starting with the μ distribution
function setμ(prim::Primitives,res::Results)
    @unpack μ, pol_func = res
    @unpack Π, A_grid, nS, nA = prim

    μ_prime = zeros(nS, nA)

    for i = 1:nS
        for j =1:nA
            A_prime = pol_func[i,j]
            for i_prime = 1:nS
                μ_prime[i_prime, j] = μ_prime[i_prime, j] + μ[i,j] * Π[i, i_prime]
            end
        end
    end
    μ_prime
end

function μ_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter
    while err>tol && n <= 1000 #begin iteration
        μ_prime = setμ(prim, res) #spit out new vectors
        err = abs.(maximum(μ_prime.-res.μ)) #reset error level
        res.μ = μ_prime #update value function
        n+=1
    end
    println("Wealth distribution converged in ", n, " iterations.")
end

## Solve the model
function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
    μ_iterate(prim, res) #iterate for μ
end
=======
#########################################
# Problem Set 2, ECON 899
# Goal Due Date: 9/20/2021
#########################################

#ReadMe: This file contains functions and variables. Run the code from the file PS2 compute.
# using Plots, Parameters, Pkg
#Keyword-enabled structure to hold model primitives

@with_kw struct Primitives
    β::Float64 = 0.9932 #discount factor
    α::Float64 = 1.5 #RRA coefficient
    e::Float64 = 1 #earnings, when employed
    u::Float64 = 0.5 #earnings, when unemployed
    nS::Int64 = 2 #number of states
    Πee::Float64 = 0.97
    Πuu::Float64 = 0.5
    Πue::Float64 = 0.03
    Πeu::Float64 = 0.5
    Π::Array{Float64, 2} = collect([Πee Πeu; Πue Πuu]) #in case I choose matrix indexing
    nA::Int64 = 1000 #number of asset grid points
    A_min::Int64 = -2 #asset lower bound
    A_max::Int64 = 5 #asset upper bound
    A_grid::Array{Float64, 1} = collect(range(A_min, length = nA, stop = A_max)) #asset grid
    a_lbar::Float64 = -2
    a_ubar::Float64 = 30
end

#Structure that holds model results
mutable struct Results
    val_func::Array{Float64, 2} #value function
    pol_func::Array{Float64, 2} #policy function
    q::Float64
    μ::Array{Float64,2}
end


#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.nS, prim.nA) #initial value function guess
    pol_func = zeros(prim.nS, prim.nA) #initial policy function guess
    q = (1 + prim.β)/2 #initial bondprice guess
    μ = zeros(prim.nS, prim.nA)
    for i = 1:prim.nS
        for j = 1:prim.nA
            μ[i,j] = (prim.A_grid[j] - prim.a_lbar)/(prim.a_ubar - prim.a_lbar)*prim.Π[i,i]
        end
    end

    res = Results(val_func, pol_func, q, μ) #initialize results struct
    prim, res #return deliverables
end

## Function for solving for the decision rule
function Bellman(prim::Primitives,res::Results)
    @unpack val_func, q = res
    @unpack β, α, e, u, Π, A_grid, nS, nA = prim
    v_next = zeros(nS, nA)

    for S_index = 1:nS
        if S_index == 1
            S = e
        else S = u
        end
        for A_index = 1:nA
            A = A_grid[A_index]
            candidate_max = -Inf
            budget = S + A
            for i = 1:nA
                c = budget - q * A_grid[i]
                if c > 0
                    val = (c^(1-α)-1)/(1 - α)
                else
                    val = -Inf
                end
                    for j = 1:nS
                        val = val + β * Π[j, S_index] * val_func[j, i]
                    end
                    if val > candidate_max
                        candidate_max = val
                        res.pol_func[S_index, A_index] = A_grid[i]
                    end
            end
            v_next[S_index, A_index] = candidate_max
        end
    end
    v_next
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter

    while err>tol && n <= 1000 #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func)) #reset error level
        res.val_func = v_next #update value function
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

## Starting with the μ distribution
function setμ(prim::Primitives,res::Results)
    @unpack μ, pol_func = res
    @unpack Π, A_grid, nS, nA = prim

    μ_prime = zeros(nS, nA)

    for i = 1:nS
        for j =1:nA
            A_prime = pol_func[i,j]
            for i_prime = 1:nS
                μ_prime[i_prime, j] = μ_prime[i_prime, j] + μ[i,j] * Π[i, i_prime]
            end
        end
    end
    μ_prime
end

function μ_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter
    while err>tol && n <= 1000 #begin iteration
        μ_prime = setμ(prim, res) #spit out new vectors
        err = abs.(maximum(μ_prime.-res.μ)) #reset error level
        res.μ = μ_prime #update value function
        n+=1
    end
    println("Wealth distribution converged in ", n, " iterations.")
end

## Finding the market clearing price
function set_price(prim::Primitives, res::Results, tol::Float64 = 1e-4)
    @unpack μ, pol_func, q = res
    xs_supply = dot(μ, pol_func)
    adj = 0.001 * q
    if abs(xs_supply) > tol && xs_supply > 0
        q_new = q - adj
        return(false)
    elseif abs(xs_supply) > tol && xs_supply < 0
        q_new = q + adj
        return(false)
    else
        println("Price is within tolerance: q = ", q)
        return(true)
    end
    res.q = q_new
end

## Solve the model
function Solve_model(prim::Primitives, res::Results, tol::Float64 = 1e-4, err::Float64 = 100.0)
    @unpack μ, pol_func = res
    converged = false
    while !converged
        V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
        μ_iterate(prim, res) #iterate for μ
        converged = set_price(prim, res)
    end
end
>>>>>>> 0621634e3580c37a0fd1a790ff6fee3604e7c60b
