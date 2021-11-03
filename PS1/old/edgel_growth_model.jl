

#keyword-enabled structure to hold model primitives
@everywhere @with_kw struct Primitives
    β::Float64 = 0.99 #discount rate
    δ::Float64 = 0.025 #depreciation rate
    θ::Float64 = 0.36 #capital share
    k_min::Float64 = 0.01 #capital lower bound
    k_max::Float64 = 75.0 #capital upper bound
    nk::Int64 = 1000 #number of capital grid points
    nz::Int64 = 2 #number of states
    k_grid::Array{Float64, 1} = collect(range(k_min, length = nk, stop = k_max)) #capital grid
    z_grid::Array{Float64, 1} = collect([1.25, 0.2])
    Π::Array{Float64, 2} = collect([0.977 0.023; 0.074 0.926])
end
#structure that holds model results
@everywhere mutable struct Results
    val_func::Array{Float64, 2} #value function
    pol_func::Array{Float64, 2} #policy function
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.nz, prim.nk) #initial value function guess
    pol_func = zeros(prim.nz, prim.nk) #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
function Bellman(prim::Primitives,res::Results)
    @unpack val_func = res #unpack value function
    @unpack k_grid, z_grid, β, δ, θ, nk, nz, Π = prim #unpack model primitives
    v_next = SharedArray{Float64}(nz, nk) #next guess of value function to fill
    k_next = SharedArray{Float64}(nz, nk)

    #choice_lower = 1 #for exploiting monotonicity of policy function
    for z_index = 1:nz
        @sync @distributed for k_index = 1:nk
            k = k_grid[k_index] #value of k
            z = z_grid[z_index] #value of z
            candidate_max = -Inf #bad candidate max
            budget = z*k^θ + (1-δ)*k #budget

            for kp_index in 1:nk #loop over possible selections of k', exploiting monotonicity of policy function
                c = budget - k_grid[kp_index] #consumption given k' selection
                if c>0 #check for positivity
                    exp_val = 0 # expected value in next period
                    for znext = 1:nz
                        exp_val = exp_val +
                                    Π[z_index, znext]*val_func[znext, kp_index]
                    end
                    val = log(c) + β*exp_val #compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        k_next[z_index, k_index] = k_grid[kp_index] #update policy function
                        #choice_lower = kp_index #update lowest possible choice
                    end
                end
            end
            v_next[z_index, k_index] = candidate_max #update value function
        end
    end
    v_next, k_next #return next guess of value function
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4,
     err::Float64 = 100.0)
    n = 0 #counter

    while err>tol #begin iteration
        v_next, k_next = Bellman(prim, res) #spit out new vectors
        err = maximum(abs.(v_next.-res.val_func)) #reset error level
        res.val_func = v_next #update value function
        res.pol_func = k_next
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

#solve the model
function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end
##############################################################################
