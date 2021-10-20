@everywhere using Parameters, DelimitedFiles, SharedArrays, LinearAlgebra

@everywhere @with_kw mutable struct  Primitives
    N_final ::Int64             = 66         # Lifespan of the agents
    J_R     ::Int64             = 46         # Retirement age
    n       ::Float64           = 0.011      # Population growth rate
    a_1     ::Float64           = 0          # Initial assets holding for newborns
    γ       ::Float64           = 0.42       # Utillity weight on consumption
    σ       ::Float64           = 2.0        # Coefficient of relative risk aversion
    α       ::Float64           = 0.36       # Capital share in production
    δ       ::Float64           = 0.06       # Capital depreciation rate
    β       ::Float64           = 0.97       # Discount factor
    N       ::Int64             = 30         # For PS4 
    θ       ::Array{Float64}    = append!([0.11], fill(0.0, N-1))

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

    # Optimal labor supply note that last argument is x = (1+r)*a-a_next
    l_opt   ::Function          = (e, w, r, a, a_next, θ) ->  (γ *(1-θ)*e*w-(1-γ)*( (1+r)*a - a_next ) ) /( (1 - θ)*w*e)

    # Production technology
    w_mkt   ::Function          = (K, L) -> (1-α)*(K^α)*(L^(-α))        # Labor first order condition
    r_mkt   ::Function          = (K, L) -> α*(K^(α-1))*(L^(1-α)) - δ   # Capital first order condition

    # Government budget constraint
    b_mkt   ::Function          = (L, w, m, θ) -> θ*w*L/m   # m is mass of retirees

    # Grids
    # Age efficiency profile
    η       ::Matrix{Float64}   = readdlm("../Data/ef.txt")
    nA      ::Int64             = 500      # Size of the asset grid
    a_min   ::Float64           = 0.0       # lower bound of the asset grid
    a_max   ::Float64           = 40.0      # upper bound of the asset grid
    a_grid  ::Array{Float64}    = collect(range(a_min, length = nA, stop = a_max))   # asset grid

end # Primitives

@everywhere mutable struct Results
    μ       ::Array{Float64, 1}             # Distibution of age cohorts
    val_fun_path ::SharedArray{Float64, 4}       # Value function for transformation
    pol_fun_path ::SharedArray{Float64, 4}       # Value function for transformation
    l_fun_path   ::SharedArray{Float64, 4}       # (effective) Labor policy function
    w_path  ::Array{Float64, 1}
    r_path  ::Array{Float64, 1}
    b_path  ::Array{Float64, 1}               
    K_path  ::Array{Float64, 1}
    L_path  ::Array{Float64, 1}
    K_path_new:: Array{Float64, 1}
    L_path_new:: Array{Float64, 1}
    F       ::Array{Float64, 4} 
end # Results

# Function that initializes the model
function Initialize()
    prim = Primitives()                   # Initialize the primitives
    d = load("../Data/Initial_Conditions.jld") # Initial Conditions from PS3
    val_fun_path = SharedArray{Float64}(prim.nA, prim.nZ,
                                        prim.N_final, prim.N)    
    pol_fun_path = SharedArray{Float64}(prim.nA, prim.nZ,
                                        prim.N_final, prim.N)    
    l_fun_path = SharedArray{Float64}(prim.nA, prim.nZ,
                                      prim.N_final, prim.N)
    K_path = collect(range(d["K_θ"], d["K"], length = prim.N))
    L_path = collect(range(d["L_θ"], d["L"], length = prim.N))
    K_path_new = zeros(prim.N) 
    L_path_new = zeros(prim.N)
    K_path_new[1] = K_path[1]
    L_path_new[1] = L_path[1]
    K_path_new[prim.N] = K_path[prim.N]
    L_path_new[prim.N] = L_path[prim.N]

    w_path = prim.w_mkt.(K_path, L_path)
    r_path = prim.r_mkt.(K_path, L_path)

    ## Population distribution across age cohorts
    μ = [1.0]
    for i in 2:prim.N_final
        push!(μ, μ[i-1]/(1.0 + prim.n))
    end
    μ = μ/sum(μ)

    b_path = prim.b_mkt.(L_path, w_path, sum(μ[prim.J_R:end]), prim.θ)

    # Initial distribution by assets, productivity level etc..
    F = SharedArray{Float64}(prim.nA, prim.nZ, prim.N_final, prim.N)    

    F[:,:,:,1] = d["Γ_0"]

    # Initialize the results
    res = Results(μ, val_fun_path, pol_fun_path, l_fun_path,
                  w_path, r_path, b_path, K_path, L_path, K_path_new,
                  L_path_new, F)

    return (prim, res)              
end


function V_ret(prim::Primitives, res::Results, t::Int64)
    @unpack nA, a_grid, N_final, J_R, util, β = prim
    @unpack b_path, r_path = res
    
    for j in N_final-1:-1:J_R
        
        choice_lower = 1

        for index_a = 1:nA

            a = a_grid[index_a]

            maxvalsofar = -Inf
            
            for index_ap = choice_lower:nA
                
                a_next = a_grid[index_ap]

                c = (1+r_path[t])*a + b_path[t] - a_next

                if c > 0
                    
                    vals = util.(c, 0) +
                        β*res.val_fun_path[index_ap, 1, j+1, t]

                    if vals > maxvalsofar
                        maxvalsofar = vals
                        res.pol_fun_path[index_a, :, j, t] .=
                            a_grid[index_ap]
                        choice_lower = index_ap
                    end

                end
            end
            res.val_fun_path[index_a, :, j, t] .= maxvalsofar
        end
    end
end



# Value function for the workers
function V_workers(prim::Primitives, res::Results, t::Int64)

    @unpack nA, nZ, z_Vals, η, N_final, J_R, util, β, θ, a_grid, Π, l_opt = prim
    @unpack r_path, w_path, b_path = res

    for j in J_R-1:-1:1 # Without progressbar for running in jupyter notebook

        # Next we iterate over the productivity levels
        for z_index in 1:nZ
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
                    l = l_opt(e, w_path[t], r_path[t], a, a_next, θ[t]) #l_grid[an_index]        # Implied labor supply in current period
                    if l < 0                    # If the labor supply is negative, we set it to zero
                        l = 0
                    elseif l > 1                # If the labor supply is greater than one, we set it to one
                        l = 1
                    end
                    c = w_path[t] * (1 - θ[t]) * e * l + (1 + r_path[t])*a - a_next # Consumption of worker (All people in this loop are working)
                    if c < 0 # If consumption is negative than this (and all future a' values) are unfeasible
                        break
                    end

                    # exp_v_next = val_fun[an_index, :, j+1] * Π[z_index , :] # Expected value of next period
                    # exp_v_next = val_fun[an_index, 1, j+1] * Π[z_index , 1] + val_fun[an_index, 2, j+1] * Π[z_index , 2] # Expected value of next period

                    # calculate expected value of next period
                    exp_v_next = 0
                    for zi = 1:nZ
                        exp_v_next += res.val_fun_path[an_index, zi, j+1, t] * Π[z_index , zi]
                    end # zi

                    v_next = util(c, l) + β * exp_v_next # next candidate to value function

                    if v_next > cand_val
                        cand_val = v_next       # Update the candidate value
                        cand_pol = a_next       # Candidate to policy function
                        cand_pol_ind = an_index # Candidate to policy function index
                        l_pol = e*l             # Candidate to labor policy function
                    end # if v_next > cand_val

                end # Next period asset choice loop

                res.val_fun_path[a_index, z_index, j, t] = cand_val         # Update the value function
                res.pol_fun_path[a_index, z_index, j, t] = cand_pol         # Update the policy function
                res.l_fun_path[a_index, z_index, j, t] = l_pol              # Update the labor policy function
                LowestChoiceInd=copy(cand_pol_ind)
            end # Current asset holdings loop
        end # Productivity loop
    end  # Age loop
end # V_workers


function ShootBackward(prim::Primitives, res::Results)
    @unpack N, util, θ, a_grid, Π, l_opt = prim
    @unpack r_path, w_path, b_path = res

    for t in N:-1:1
        
        last_period_value = util.( a_grid .*(1  + r_path[t]) .+ b_path[t], 0)
    
        res.val_fun_path[: ,: , end, t] = hcat(last_period_value,
                                               last_period_value)
        V_ret(prim, res, t);
        
        V_workers(prim, res, t);

    end
    return res
end


function CalculatedDist(prim::Primitives, res::Results, t::Int64)

    @unpack N_final, n, p_L, p_H, nZ, nA, Π, a_grid, N = prim
    
    # Finding the steady state distribution
    for j in 2:N_final
        for z_ind in 1:nZ
            for a_ind in 1:nA
                a_next_ind = argmin(abs.(res.pol_fun_path[a_ind, z_ind, j-1, t-1].-a_grid))
                
                for zi = 1:nZ
                    res.F[a_next_ind, zi, j, t] +=
                        res.F[a_ind, z_ind, j-1, t-1] * Π[z_ind, zi] * (res.μ[j]/res.μ[j-1])
                end # zi
            end
        end # z_ind
        res.F[1, 1, 1, t] = res.μ[1] * p_H ## New born each period
        res.F[1, 2, 1, t] = res.μ[1] * p_L ## New born each period
    end # j loop
end

function UpdatePath(prim::Primitives, res::Results, t::Int64)

    @unpack a_grid = prim
    
    res.K_path_new[t] = sum(res.F[:,:,:,t].* a_grid)
    res.L_path_new[t] = sum(res.F[:,:,:,t].* res.l_fun_path[:,:,:,t])

end

function ShootForward(prim::Primitives, res::Results)

    for t = 2:prim.N-1
        CalculatedDist(prim, res, t)
        UpdatePath(prim, res, t)
    end
end

function Convergence()
                   
    prim, res = Initialize();

    @unpack r_mkt, w_mkt, b_mkt, J_R, θ = prim
    err = 100;
    λ = 0.95;
    for n = 1:10
        println("$n iterations; err = $err")

        res.r_path = r_mkt.(res.K_path, res.L_path)
        res.w_path = w_mkt.(res.K_path, res.L_path)
        res.b_path = b_mkt.(res.L_path, res.w_path,
                            sum(res.μ[J_R:end]), θ)
        ShootBackward(prim, res);
        ShootForward(prim, res);

        err = maximum(abs.(res.K_path_new - res.K_path) +
                      abs.(res.L_path_new - res.L_path))

        if err > 1e-2 
            res.K_path = λ.*res.K_path_new + (1 - λ).*res.K_path
            res.L_path = λ.*res.L_path_new + (1 - λ).*res.L_path
        end
    end
    return prim, res
end



