
@everywhere using Parameters, DelimitedFiles, ProgressBars, SharedArrays, LinearAlgebra

# Define the primitives of the model
@everywhere @with_kw mutable struct  Primitives
    N_final ::Int64             = 66         # Lifespan of the agents
    J_R     ::Int64             = 46         # Retirement age
    n       ::Float64           = 0.011      # Population growth rate
    a_1     ::Float64           = 0          # Initial assets holding for newborns
    θ       ::Float64           = 0.11       # Labor income tax rate
    γ       ::Float64           = 0.42       # Utillity weight on consumption
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

    util  ::Function          = (c, l) -> ( c > 0 ) ? ((c^γ * (1 - l)^γ)^(1-σ))/(1-σ) : -Inf

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
    η       ::Matrix{Float64}   = readdlm("./PS3/Data/ef.txt")
    nA      ::Int64             = 1000      # Size of the asset grid
    a_min   ::Float64           = 0.0       # lower bound of the asset grid
    a_max   ::Float64           = 75.0      # upper bound of the asset grid
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

    # ! This is a experiment, maybe it is useful to also save
    # ! the indices of the optimal policy function
    pol_fun_ind ::Array{Int64,3}            # Policy function indices
    F       ::Array{Float64,3}              # Distribution of agents over asset holdings
end # Results

# Function that initializes the model
function Initialize()
    prim = Primitives()                             # Initialize the primitives
    w = 1.05                                        # Wage guess
    r = 0.05                                        # Interest rate guess
    b = 0.2                                         # Benefits guess
    K = 4                                           # inital capital guess
    L = 0.9                                         # initial labor guess
    val_fun = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final)    # Initialize the value function
    pol_fun = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final)    # Initialize the policy function
    pol_fun_ind = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final)# Initialize the policy function indices
    l_fun = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final)

    # Before the model starts, we can set the initial value function at the end stage
    # We set the last age group to have a value function consuming all the assets and
    # with a labor supply 0 (i.e. no labor) and recieving a benefit of b
    last_period_value = prim.util.( prim.a_grid .* (1 + r) .+ b, 0 )
    val_fun[: ,: , end] = hcat(last_period_value, last_period_value)

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
    res = Results(w, r, b, K, L, μ, val_fun, pol_fun, l_fun, pol_fun_ind, F)

    return (prim, res)                              # Return the primitives and results
end

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
            res.pol_fun_ind[a_index, :, j] .= pol_ind
            res.pol_fun[a_index, :, j] .= a_grid[pol_ind]
            res.val_fun[a_index, :, j] .= val_max
            res.l_fun[a_index, :, j] .= 0
        end # for a_index
    end # for j

end # V_ret

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
                res.pol_fun_ind[a_index, z_index, j] = cand_pol_ind # Update the policy function index
                res.l_fun[a_index, z_index, j] = l_pol              # Update the labor policy function
                LowestChoiceInd=copy(cand_pol_ind)
            end # Current asset holdings loop
        end # Productivity loop
    end  # Age loop
    #=
    res.val_fun = val_fun
    res.pol_fun = pol_fun
    res.pol_fun_ind = pol_fun_ind
    res.l_fun = l_fun
    =#
end # V_workers

# If we want to speed up the code we can use Fortran
# the following function is a wrapper for the Fortran code

function V_Fortran(prim::Primitives, res::Results)
    @unpack r, w, b =res
    # PS3/FortranCode/conesa_kueger.f90
    # Compile Fortran code
    path = "./PS3/FortranCode/"
    run(`gfortran -fopenmp -O2 -o $(path)V_Fortran $(path)conesa_kueger.f90`)
    # run(`./T_op $q $n_iter`)
    run(`$(path)V_Fortran`)

    results_raw =  readdlm("$(path)results.csv");

    val_fun = zeros(prim.nA, prim.nZ, prim.N_final)    # Initialize the value function
    pol_fun = zeros(prim.nA, prim.nZ, prim.N_final)    # Initialize the policy function
    pol_fun_ind = zeros(prim.nA, prim.nZ, prim.N_final + 1)    # Initialize the policy function index
    consumption = zeros(prim.nA, prim.nZ, prim.N_final)    # Initialize the consumption function
    l_fun = zeros(prim.nA, prim.nZ, prim.N_final)    # Initialize the labor policy function
    for j in 1:prim.N_final
        range_a = (j-1) * 2*prim.nA + 1 : j * 2*prim.nA |> collect
        val_fun[:,:,j] = hcat(results_raw[range_a[1:prim.nA],end], results_raw[range_a[prim.nA+1:end],end])
        pol_fun_ind[:,:,j] = hcat(results_raw[range_a[1:prim.nA],end-1], results_raw[range_a[prim.nA+1:end],end-1])
        pol_fun[:,:,j] = hcat(results_raw[range_a[1:prim.nA],end-2], results_raw[range_a[prim.nA+1:end],end-2])
        consumption[:,:,j]   = hcat(results_raw[range_a[1:prim.nA],end-3], results_raw[range_a[prim.nA+1:end],end-3])
        l_fun[:,:,j] = hcat(results_raw[range_a[1:prim.nA],end-4], results_raw[range_a[prim.nA+1:end],end-4])
    end
    A_grid_fortran = results_raw[1:prim.nA,end-5]
    res.val_fun = val_fun
    res.pol_fun = pol_fun
    res.pol_fun_ind = pol_fun_ind
    res.l_fun = l_fun
    return A_grid_fortran, consumption

end # run_Fortran()


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

    prim, res= Initialize()

    # make relevant policy experiment changes
    if ~ss 
        prim.θ = 0
    end

    if ~i_risk 
        prim.z_H, prim.z_L, prim.z_Vals = 0.5, 0.5, [0.5, 0.5]
    end

    if exog_l
        prim.γ = 1
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
            A_grid_fortran, consumption = V_Fortran(prim,res)
        else
            V_ret(prim, res);
            V_workers(prim, res);
        end
        SteadyStateDist(prim, res);

        # calculate aggregate capital and labor
        K = sum(res.F[:, :, :] .* a_grid)
        L = sum(res.F[:, :, :] .* res.l_fun) # Labor supply grid

        # calculate error
        err = maximum(abs.([res.K, res.L] - [K, L]))

        if (err > tol*10)
            # Leave λ at the default
        elseif (err > tol*5) & (λ <= 0.85)
            λ = 0.85
        elseif (err > tol*1.01) & (λ <= 0.90)
            λ = 0.90
        elseif λ <= 0.95
            λ = 0.95
        end

        # update guess
        res.K = (1-λ)*K + λ*res.K
        res.L = (1-λ)*L + λ*res.L

        n+=1

        println("$n iterations; err = $err, K = ", round(res.K, digits = 4), ", L = ",
        round(res.L, digits = 4), ", λ = $λ")

    end # while err > tol
    return prim, res
end # MarketClearing

# Function to calculate compensating variation
function Lambda(prim::Primitives, res::Results, W::Float64)
    
    # unpack necessary variables
    @unpack F, val_fun = res
    @unpack α, β = prim

    # calculate and return compensating variation
    λ = ( (W + (1/((1-α)*(1-β)))) ./ (val_fun .+ (1/((1-α)*(1-β)))) ).^(1/(1-α)) .- 1

    return F.*λ

end
