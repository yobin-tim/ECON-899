
@everywhere using Parameters, DelimitedFiles, ProgressBars, SharedArrays, LinearAlgebra, NaNMath

# Define the primitives of the model
@everywhere @with_kw mutable struct  Primitives
    N_final ::Int64             = 66         # Lifespan of the agents
    J_R     ::Int64             = 46         # Retirement age
    n       ::Float64           = 0.011      # Population growth rate
    a_1     ::Float64           = 0          # Initial assets holding for newborns
    θ       ::Float64                        # Labor income tax rate
    γ       ::Float64                        # Utillity weight on consumption
    σ       ::Float64           = 2.0        # Coefficient of relative risk aversion
    α       ::Float64           = 0.36       # Capital share in production
    δ       ::Float64           = 0.06       # Capital depreciation rate
    β       ::Float64           = 0.97       # Discount factor

    # Parameters regarding stochastic processes
    z_H     ::Float64           = 3.0        # Idiosyncratic productivity High
    z_L     ::Float64           = 0.5        # Idiosyncratic productivity Low
    z_Vals  ::Array{Float64}    = [z_H, z_L] # Vector of idiosyncratic productivity values
    nZ      ::Int64             = 2          # Number of idiosynctatic productivity levels
    p_H     ::Float64           = 0.2037     # Probability of z_H at birth
    p_L     ::Float64           = 0.7963     # Probability of z_L at birth

    # Markov transition matrix for z
    Π       ::Array{Float64,2}  = [0.9261 1-0.9261; 1- 0.9811  0.9811]

    # Functions

    util  ::Function          = (c, l) -> ( c > 0 ) ? ((c^γ * (1 - l)^(1-γ))^(1-σ))/(1-σ) : -Inf

    # Utility of a retiree
    # Todo: Remove the next 3 lines if everything is working
    # * Note im only using the utility of a worker and setign l = 0 to obtain the utility of a retiree
    # util_R  ::Function          = (c) -> c^(γ*(1-σ))/(1-σ)

    # Optimal labor supply note that last argument is x = (1+r)*a-a_next
    l_opt   ::Function          = (e, w, r, a, a_next) ->  (γ *(1-θ)*e*w-(1-γ)*( (1+r)*a - a_next ) ) /( (1 - θ)*w*e)

    # Production technology
    w_mkt   ::Function          = (K, L) -> (1-α)*(K^α)*(L^(-α))        # Labor first order condition
    r_mkt   ::Function          = (K, L) -> α*(K^(α-1))*(L^(1-α)) - δ   # Capital first order condition

    # Government budget constraint
    b_mkt   ::Function          = (L, w, m) -> θ*w*L/m   # m is mass of retirees

    # Grids
    # Age efficiency profile
    η       ::Matrix{Float64}   = readdlm("../Data/ef.txt")
    nA      ::Int64             = 500      # Size of the asset grid
    a_min   ::Float64           = 0.0       # lower bound of the asset grid
    a_max   ::Float64           = 40.0      # upper bound of the asset grid
    a_grid  ::Array{Float64}    = collect(range(a_min, length = nA, stop = a_max))   # asset grid

end # Primitives

# Structure mutable parameters and results of the model
@everywhere mutable struct Results
    w       ::Float64                       # Wage
    r       ::Float64                       # Interest rate
    b       ::Float64                       # Benefits
    K       ::Float64                       # aggregate capital
    L       ::Float64                       # aggregate labor
    μ       ::Array{Float64, 1}             # Distibution of age cohorts
    val_fun ::SharedArray{Float64, 3}             # Value function
    pol_fun ::SharedArray{Float64, 3}             # Policy function
    l_fun   ::SharedArray{Float64, 3}             # (effective) Labor policy function

    F       ::Array{Float64,3}              # Distribution of agents over asset holdings
end # Results

# Function that initializes the model
function Initialize(; θ = 0.11, γ = 0.42)
    prim = Primitives(θ = θ, γ = γ)                 # Initialize the primitives
    w = 1.05                                        # Wage guess
    r = 0.05                                        # Interest rate guess
    b = 0.2                                         # Benefits guess
    K = 4                                           # inital capital guess
    L = 0.9                                         # initial labor guess
    val_fun = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final)    # Initialize the value function
    pol_fun = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final)    # Initialize the policy function
    l_fun = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final)

    # Before the model starts, we can set the initial value function at the end stage
    # We set the last age group to have a value function consuming all the assets and
    # with a labor supply 0 (i.e. no labor) and recieving a benefit of b

    # Calculate population distribution across age cohorts
    μ = [1.0]
    for i in 2:prim.N_final
        push!(μ, μ[i-1]/(1.0 + prim.n))
    end
    μ = μ/sum(μ)

    # Finally we initialize the distribution of the agents
    F = zeros(prim.nA, prim.nZ, prim.N_final)
    F[1, 1, 1] = μ[1] * prim.p_H
    F[1, 2, 1] = μ[1] * prim.p_L

    # Initialize the results
    res = Results(w, r, b, K, L, μ, val_fun, pol_fun, l_fun, F)

    return (prim, res)                              # Return the primitives and results
end
#=
# Value funtion for the retirees
function V_ret(prim::Primitives, res::Results)

    # unpack the primitives and the results
    @unpack nA, a_grid, N_final, J_R, util = prim
    @unpack b, r = res

    # We obtain for every age group and asset holdings level the value function using backward induction
    for j in N_final-1:-1:J_R
        @sync @distributed for a_index in 1:nA
            a = a_grid[a_index]
            vals = util.(((1+r)*a + b).- a_grid, 0) .+ res.val_fun[:, 1, j+1]
            pol_ind = argmax(vals)
            val_max = vals[pol_ind]
            res.pol_fun[a_index, :, j] .= a_grid[pol_ind]
            res.val_fun[a_index, :, j] .= val_max
            res.l_fun[a_index, :, j] .= 0
        end # for a_index
    end # for j

end # V_ret
=#
function V_last(prim::Primitives, res::Results)
    last_period_value = prim.util.( prim.a_grid .* (1 + res.r) .+ res.b, 0 )
    res.val_fun[: ,: , end] = hcat(last_period_value, last_period_value)
end

function V_ret(prim::Primitives, res::Results)
    @unpack nA, a_grid, N_final, J_R, util, β = prim
    @unpack b, r = res

    for j in N_final-1:-1:J_R
        
        choice_lower = 1

        for index_a = 1:nA

            a = a_grid[index_a]

            maxvalsofar = -Inf
            
            for index_ap = choice_lower:nA
                
                a_next = a_grid[index_ap]

                c = (1+r)*a+b - a_next

                if c > 0
                    
                    vals = util.(c, 0) +
                        β*res.val_fun[index_ap, 1, j+1]

                    if vals > maxvalsofar
                        maxvalsofar = vals
                        res.pol_fun[index_a, :, j] .=
                            a_grid[index_ap]
                        choice_lower = index_ap
                    end

                end
            end
            res.val_fun[index_a, :, j] .= maxvalsofar
        end
    end
end



# Value function for the workers
function V_workers(prim::Primitives, res::Results)

    # Unopack the primitives
    @unpack nA, nZ, z_Vals, η, N_final, J_R, util, β, θ, a_grid, Π, l_opt = prim
    @unpack r, w, b = res

    # First  we iterate over the age groups
    # for j in ProgressBar(J_R-1:-1:1) # Progressbar for running in console

    for j in J_R-1:-1:1 # Without progressbar for running in jupyter notebook

        # Next we iterate over the productivity levels
        @sync @distributed for z_index in 1:nZ
            z = z_Vals[z_index] # Current idiosyncratic productivity level
            #println("Solving for productivity type $z")
            e = z * η[j] # Worker productivity level (only for working age)
            LowestChoiceInd=1 #Exploiting monotonicity in the policy function
            # Next we iterate over the asset grid
            for a_index in 1:nA
                a = a_grid[a_index] # Current asset level
                cand_val = -Inf     # Initialize the candidate value
                cand_pol = 0        # Initialize the candidate policy
                cand_pol_ind = 1    # Initialize the candidate policy index
                l_pol = 0           # Initialize the labor policy

                # Next we iterate over the possible choices of the next period's asset
                #l_grid = l_opt.(e, w, r, a, a_grid) # Labor supply grid
                # if j == 20 && z_index == 2
                #     print("\n a = $a a_next reached:")
                # end
                for an_index in LowestChoiceInd:nA
                    a_next = a_grid[an_index]   # Next period's asset level
                    l = l_opt(e, w, r, a, a_next) #l_grid[an_index]        # Implied labor supply in current period
                    if l < 0                    # If the labor supply is negative, we set it to zero
                        l = 0
                    elseif l > 1                # If the labor supply is greater than one, we set it to one
                        l = 1
                    end
                    c = w * (1 - θ) * e * l + (1 + r)*a - a_next # Consumption of worker (All people in this loop are working)
                    if c < 0 # If consumption is negative than this (and all future a' values) are unfeasible
                        break
                    end

                    # exp_v_next = val_fun[an_index, :, j+1] * Π[z_index , :] # Expected value of next period
                    # exp_v_next = val_fun[an_index, 1, j+1] * Π[z_index , 1] + val_fun[an_index, 2, j+1] * Π[z_index , 2] # Expected value of next period

                    # calculate expected value of next period
                    exp_v_next = 0
                    for zi = 1:nZ
                        exp_v_next += res.val_fun[an_index, zi, j+1] * Π[z_index , zi]
                    end # zi

                    v_next = util(c, l) + β * exp_v_next # next candidate to value function

                    if v_next > cand_val
                        cand_val = v_next       # Update the candidate value
                        cand_pol = a_next       # Candidate to policy function
                        cand_pol_ind = an_index # Candidate to policy function index
                        l_pol = e*l             # Candidate to labor policy function
                    end # if v_next > cand_val

                end # Next period asset choice loop

                res.val_fun[a_index, z_index, j] = cand_val         # Update the value function
                res.pol_fun[a_index, z_index, j] = cand_pol         # Update the policy function
                res.l_fun[a_index, z_index, j] = l_pol              # Update the labor policy function
                LowestChoiceInd=copy(cand_pol_ind)
            end # Current asset holdings loop
        end # Productivity loop
    end  # Age loop
end # V_workers

# Function to obtain the steady state distribution
function SteadyStateDist(prim::Primitives, res::Results)
    # Initialize the steady state distribution
    res.F[:,:,2:end] .= zeros(prim.nA, prim.nZ)
    # Unpack the primitives
    @unpack N_final, n, p_L, p_H, nZ, nA, Π, a_grid = prim

    # Finding relative size of each age cohort

    # Finding the steady state distribution
    for j in 2:N_final
        for z_ind in 1:nZ
            for a_ind in 1:nA
                a_next_ind = argmin(abs.(res.pol_fun[a_ind, z_ind, j-1].-a_grid))
                #= This should not ever happen
                if a_next_ind == 0 # Level not reached
                    continue
                end
                =#

                for zi = 1:nZ
                    res.F[a_next_ind, zi, j] += res.F[a_ind, z_ind, j-1] * Π[z_ind, zi] * (res.μ[j]/res.μ[j-1])
                end # zi

            end
        end # z_ind
    end # j loop

end # SteadyStateDist

# Function to solve for market prices
function MarketClearing(; ss::Bool=true, i_risk::Bool=true, exog_l::Bool=false,
    use_Fortran::Bool=false, λ::Float64=0.7, tol::Float64=1e-2, err::Float64=100.0)

    # initialize struct according to policies
    if ~ss & exog_l
        prim, res = Initialize(θ = 0, γ = 1)
        prim.l_opt = (e, w, r, a, a_next) ->  1
        res.L = sum(res.μ[1:(prim.J_R-1)])
    elseif ~ss
        prim, res = Initialize(θ = 0)
    elseif exog_l 
        prim, res = Initialize(γ = 1)
        prim.l_opt = (e, w, r, a, a_next) ->  1
        res.L = sum(res.μ[1:(prim.J_R-1)])
    else
        prim, res = Initialize()
    end

    if ~i_risk 
        prim.z_H, prim.z_L, prim.z_Vals = 0.5, 0.5, [0.5, 0.5]
    end

    # unpack relevant variables and functions
    @unpack w_mkt, r_mkt, b_mkt, J_R, a_grid = prim

    n = 0 # loop counter

    # iteratively solve the model until excess savings converge to zero
    while err > tol

        # calculate prices and payments at current K, L, and F
        res.r = r_mkt(res.K, res.L)
        res.w = w_mkt(res.K, res.L)
        res.b = b_mkt(res.L, res.w, sum(res.μ[J_R:end]))

        # solve model with current model and payments
        if use_Fortran
            K, L = V_Fortran(prim,res)
        else
            V_last(prim, res);
            V_ret(prim, res);
            V_workers(prim, res);
            SteadyStateDist(prim, res);
    
            # calculate aggregate capital and labor
            K = sum(res.F[:, :, :] .* a_grid)
            L = sum(res.F[:, :, :] .* res.l_fun) # Labor supply grid
        end

        # calculate error
        err = maximum(abs.([res.K, res.L] - [K, L]))

        if (err > tol*10)
            # Leave λ at the default
        elseif (err > tol*5) & (λ <= 0.85)
            λ = 0.85
        elseif (err > tol*1.7) & (λ <= 0.90)
            λ = 0.90
        elseif λ <= 0.975
            λ = 0.975
        end

        # update guess
        res.K = (1-λ)*K + λ*res.K
        res.L = (1-λ)*L + λ*res.L

        n+=1

        if  n % 5 == 0
            println("$n iterations; err = $err, K = ", round(res.K, digits = 4), ", L = ",
            round(res.L, digits = 4), ", λ = $λ")
        end

    end # while err > tol
    return prim, res
end # MarketClearing

# Function to calculate compensating variation
function Lambda(prim::Primitives, res::Results, W::SharedArray{Float64, 3})
    
    # unpack necessary variables
    @unpack F, val_fun = res
    @unpack α, β, γ, σ = prim

    # calculate and return compensating variation
    # λ = (W ./ val_fun).^(1/(γ*(1-σ))) .- 1

    λ = (val_fun ./ W).^(1/(γ*(1-σ))) .- 1
    # Is W denominator?

    return NaNMath.sum(F.*λ)

end

# To anwer "who benefits" in Exercise 3.  
function Lambda2(prim::Primitives, res::Results, W::SharedArray{Float64, 3})
    
    # unpack necessary variables
    @unpack F, val_fun = res
    @unpack α, β, γ, σ = prim

    # calculate and return compensating variation
    λ = (val_fun ./ W).^(1/(γ*(1-σ))) .- 1
    a = [1:1:66;]
    b = dropdims(sum(sum(F.*λ, dims = 1), dims = 2), dims =1)'
    c = [a, b]
    return c     
end
