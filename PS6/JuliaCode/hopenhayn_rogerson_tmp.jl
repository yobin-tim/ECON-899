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
    p_max       ::Float64 = 3.0
    # nP          ::Int64   = 10
    # p_grid      ::Array{Float64}   = range(p_min, stop = p_max, length = nP)

    # Optimal decision rules
    n_optim     ::Function         = ( s , p ) -> (θ * p * s) ^ (1/(1 - θ))

    Π           ::Function         = ( s , p,  n ) -> (n > 0) ? p*s*( n )^θ - n  - p * c_f : -p * c_f

    # Limits for the mass of entrants
    M_min       ::Float64 = 1.0
    M_max       ::Float64 = 10.0


end

# Structure that stores results
mutable struct Results
    W_val   ::Array{Float64}         # Firm value given state variables
    n_opt   ::Array{Float64}         # Optimal labor demand for each possible state
    x_opt   ::Array{Float64}           # Optimal firm decicion for each state
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
    p = (prim.p_max + prim.p_min)/2
    μ = ones(prim.nS) / prim.nS # Uniform distribution is the initial guess
    M = 5.0

    res = Results(W_val, n_opt, x_opt, p, μ, M)

    return prim, res

end

function 
