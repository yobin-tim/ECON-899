using Parameters, DelimitedFiles, ProgressBars
# Define the primitives of the model
@with_kw mutable struct  Primitives
    N_final ::Int64             = 66         # Lifespan of the agents
    J_R     ::Int64             = 46         # Retirement age
    n       ::Float64           = 0.011      # Population growth rate
    a_1     ::Float64           = 0          # Initial assets holding for newborns
    θ       ::Float64           = 0.11       # Labor income tax rate
    γ       ::Float64           = 0.42       # Utillity weight on consumption
    σ       ::Float64           = 2.0        # Coefcient of relative risk aversion
    α       ::Float64           = 0.36       # Capital share in production
    δ       ::Float64           = 0.06       # Capital depreciation rate
    β       ::Float64           = 0.97       # Discount factor
    
    # Parameters regarding stochastic processes
    z_H     ::Float64           = 1.0        # Idiosyncratic productivity High
    z_L     ::Float64           = 0.5        # Idiosyncratic productivity Low
    z_Vals  ::Array{Float64}    = [z_H, z_L] # Vector of idiosyncratic productivity values
    nZ      ::Int64             = 2          # Number of idiosynctatic productivity levels
    p_H     ::Float64           = 0.2037     # Probability of z_H at birth
    p_L     ::Float64           = 0.7963     # Probability of z_L at birth
    # Markov transition matrix for z
    Π       ::Array{Float64,2}  = [0.9261 1-0.9261; 1- 0.9811  0.9811] 

    # Functions
    # Utility of a worker
    # * Previously calle util_w 
    # Todo: Change name to util_w if there is a problem
    util  ::Function          = (c, l) -> ( c > 0 ) ? (c^γ * (1 - l)^γ)^(1-σ)/(1-σ) : -Inf
    # Utility of a retiree
    # Todo: Remove the next 3 lines if everything is working
    # * Note im only using the utility of a worker and setign l = 0 to obtain the utility of a retiree
    # util_R  ::Function          = (c) -> c^(γ*(1-σ))/(1-σ)

    # Optimal labor supply note that last argument is x = (1+r)*a-a_next
    l_opt   ::Function          = (e, w, r, a, a_next) ->  (γ *(1-θ)*e*w-(1-γ)*( (1+r)*a - a_next ) ) /( (1 - θ)*w*e) 
    
    # Grids
    # Age efficiency profile
    η       ::Matrix{Float64}   = readdlm("./PS3/Data/ef.txt") 
    nA      ::Int64             = 1000      # Size of the asset grid
    a_min   ::Float64           = 0.0       # lower bound of the asset grid
    a_max   ::Float64           = 100.0      # upper bound of the asset grid
    a_grid  ::Array{Float64}    = collect(range(a_min, length = nA, stop = a_max))   # asset grid
end # Primitives

# Structure mutable parameters and results of the model 
mutable struct Results
    w       ::Float64                       # Wage        
    r       ::Float64                       # Interest rate
    b       ::Float64                       # Benefits
    l_grid ::Array{Float64, 4}              # Labor supply grid
    val_fun ::Array{Float64,3}              # Value function
    pol_fun ::Array{Float64,3}              # Policy function
    # ! This is a experiment, maybe it is usefull to also save 
    # ! the indices of the optimal policy function
    pol_fun_ind ::Array{Int64,3}            # Policy function indices
end # Results

# Function that initializes the model
function Initialize()
    prim = Primitives()                             # Initialize the primitives
    w = 1.05                                        # Wage
    r = 0.05                                        # Interest rate
    b = 0.2                                         # Benefits
    val_fun = zeros(prim.nA, prim.nZ, prim.N_final)    # Initialize the value function
    pol_fun = zeros(prim.nA, prim.nZ, prim.N_final)    # Initialize the policy function
    pol_fun_ind = zeros(prim.nA, prim.nZ, prim.N_final)# Initialize the policy function indices
    
    # ! This is not showing much improvement in speed
    # We can pre compute the labor supply grid
    l_grid = zeros(prim.nZ, prim.N_final, prim.nA, prim.nA)    # labor grid
    # for z in 1:prim.nZ, age in 1:prim.N_final, a in 1:prim.nA, an in 1:prim.nA
    #     if age < prim.J_R
    #         x = (1+r)*prim.a_grid[a] - prim.a_grid[an]
    #         e = prim.η[age] * prim.z_Vals[z]
    #         l_grid[z, age, a, an] = prim.l_opt(e, w, x)
    #     else
    #         l_grid[z, age, a, an] = 0
    #     end # if
    # end # for
    
    # Before the model starts, we can set the initial value function at the end stage
    # We set the last age group to have a value function consuming all the assets and
    # with a labor supply 0 (i.e. no labor) and recieving a benefit of b
    last_period_value = prim.util.( prim.a_grid .* (1 + r) .+ b, 0 ) 
    val_fun[: ,: , end] = hcat(last_period_value, last_period_value) 
    
    # Initialize the results
    res = Results(w, r, b, l_grid, val_fun, pol_fun, pol_fun_ind)        
    
    return (prim, res)                              # Return the primitives and results
end

# Value funtion for the retirees
function V_ret(prim::Primitives, res::Results)
    # unpack the primitives and the results
    @unpack nA, a_grid, N_final, J_R, util = prim
    @unpack b, r = res
    # We obtain for every age group and asset holdigs level the value function using backward induction
    for j in N_final-1:-1:J_R
        for a_index in 1:nA
            a = a_grid[a_index]
            vals = util.(((1+r)*a + b ).- a_grid, 0) .- res.val_fun[:, 1, j+1]
            pol_ind = argmax(vals)
            val_max = vals[pol_ind]
            res.pol_fun_ind[a_index, :, j] .= pol_ind
            res.pol_fun[a_index, :, j] .= a_grid[pol_ind]
            res.val_fun[a_index, :, j] .= val_max
        end # for a_index 
    end # for j
end # V_ret

# Value function for the workers
function V_workers(prim::Primitives, res::Results)
    # Unopack the primitives
    @unpack nA, nZ, z_Vals, η, N_final, J_R, util, β, θ, a_grid, Π, l_opt = prim
    @unpack r, w, b, l_grid, val_fun, pol_fun, pol_fun_ind = res

    # First we iterate over the productivity levels
    for z_index in 1:nZ
        z = z_Vals[z_index] # Current idiosyncratic productivity level
        println("Solving for productivity type $z")
        # Next we iterate over the age groups
        for j in ProgressBar(J_R-1:-1:1) # Progressbar for runing in console
        # for j in N_final-1:-1:1 # Without progressbar for runing in jupyter notebook
            e = ( j < J_R ) ? z * η[j] : 0 # Worker productivity level (only for working age)
            # Next we iterate over the asset grid
            for a_index in 1:nA
                a = a_grid[a_index] # Current asset level
                cand_val = -Inf # Initialize the candidate value
                cand_pol = 0 # Initialize the candidate policy
                cand_pol_ind = 0 # Initialize the candidate policy index
                # Next we iterate over the possible choices of the next period's asset
                l_grid = l_opt.(e, w, r, a, a_grid) # Labor supply grid
                for an_index in 1:nA
                    a_next = a_grid[an_index] # Next period's asset level
                    l = l_grid[an_index] # Current labor supply
                    if ( j < J_R ) # If the agent is working
                        c = w * (1 - θ) * e * l + (1 + r)a - a_next # Consumption
                    else # If the agent is not working
                        c = (1 + r) * a - a_next # Consumption of retiree
                    end
                    if c < 0 # If the consumption is negative, stop the exploration
                        break
                    end
                    # exp_v_next = val_fun[an_index, :, j+1] * Π[z_index , :] # Expected value of next period
                    exp_v_next = val_fun[an_index, 1, j+1] * Π[z_index , 1] + val_fun[an_index, 2, j+1] * Π[z_index , 2] # Expected value of next period
                    v_next = util(c, l) + β * exp_v_next # next candidate to value function
                    if v_next > cand_val
                        cand_val = v_next # Update the candidate value
                        cand_pol = a_next # Candidate to policy function
                        cand_pol_ind = an_index # Candidate to policy function index
                    end # if
                end # Next period asset choice loop
                val_fun[a_index, z_index, j] = cand_val # Update the value function
                pol_fun[a_index, z_index, j] = cand_pol # Update the policy function
                pol_fun_ind[a_index, z_index, j] = cand_pol_ind # Update the policy function index
            end # Current asset holdings loop
        end # Age loop
    end # Productivity loop
    res.val_fun = val_fun
    res.pol_fun = pol_fun
    res.pol_fun_ind = pol_fun_ind
end # V_workers

# Function to obtain the steady state distribution
function SteadyStateDist(prim::Primitives, res::Results)
    # Unpack the primitives
    @unpack N_final, n = prim
    # Finding relative size of each age cohort
    μ = [1.0]
    for i in 2:N_final
        push!(μ, μ[i-1]/(1.0 + n))
    end

    μ/sum(μ)
end

# SteadyStateDist(prim, res)