using Parameters

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
    p_max       ::Float64 = 15.0
    # nP          ::Int64   = 10
    # p_grid      ::Array{Float64}   = range(p_min, stop = p_max, length = nP)

    # Optimal decision rules
    n_optim     ::Function         = ( s , p ) -> (θ * p * s) ^ (1/(1 - θ))

    Π           ::Function         = ( s , p ) -> p*s*( n_optim( s , p) )^θ - n_optim( s , p) - p * c_f

end

# Structure that stores results
mutable struct Results
    W_val   ::Array{Float64}         # Firm value given state variables
    n_opt   ::Array{Float64}         # Optimal labor demand for each possible state
    x_opt   ::Array{Int64}           # Optimal firm decicion for each state
    p       ::Float64                # Price that clears the market
end

# Initialize model
function Initialize()
    prim  = Primitives()
    
    W_val = zeros(prim.nS)
    n_opt = zeros(prim.nS)
    x_opt = zeros(prim.nS)
    p = 7.0

    res = Results(W_val, n_opt, x_opt, p)

    return prim, res

end

# Bellman operator for W
function W(p, Π, nS, trans_mat)
        
    # Iterate over all possible states
    for s_i ∈ 1:nS

        prof = Π(s_i, p)

        # Calculate expected continuation value

        exp_cont_value = trans_mat[s_i, :]' * res.W_val

        x = ( exp_cont_value > 0 ) ? 1 : 0
        
        res.W_val[s_i] = prof + (1 - x) * exp_cont_value
        res.x_opt[s_i] = x

    end

end # W

#Value function iteration for W operator
function TW_iterate(p, Π, nS, trans_mat; tol::Float64 = 1e-4)
    n = 0 #counter
    err = 100.0 #initialize error
    while  (err > tol) & (n < 4000)#begin iteration
        W_val_old = copy(res.W_val)
        W(p, Π, nS, trans_mat)
        err = abs.( maximum( W_val_old - res.W_val ) ) #reset error level
        n+=1
        if n % 100 == 0
            println("Iter =", n , " Error = ", err)
        end
    end
end # TW_iterate

function market_clearing(prim, res; tol::Float64 = 1e-3)
    @unpack Π, nS, trans_mat, ν, c_e, p_min, p_max = prim

    θ = 0.99
    n = 0
    
    while n < 1000
        TW_iterate(res.p, Π, nS, trans_mat)
        # Calculate EC
        EC = sum(res.W_val .* ν)/res.p - c_e
        # println("p = ", res.p," EC = ", EC, " tol = ", tol)
        println(n, " iterations; EC = ", EC, ", p = ", res.p, ", θ = ", θ)
        if abs( EC ) < tol
            println("Market Cleared")
            break
        end
        if EC > 0       # adjust price toward bounds according to tuning parameter
            res.p = θ*res.p + (1-θ)*p_min
        else
            res.p = θ*res.p + (1-θ)*p_max
        end
        n+=1


        # adjust tuning parameter based on EC
        if abs(EC) < tol * 1000
            θ = 0.999
        elseif abs(EC) > tol * 10000
            θ = 0.95
        else
            θ = 0.99
        end
        
        n += 1
        
    end

end

