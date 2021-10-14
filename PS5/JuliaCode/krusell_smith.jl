using LinearAlgebra, Parameters
@doc """
    The following function recieves the followig input:

        - d_z:Array{Float64,1} the duration of states [d_g, d_b]
        - d_unemp:Array{Float64,1} the duration of unemployment in each state [d_unemp_g, d_unemp_b]
        - u:Array{Float64,1} the fraction of people unemployed inn each state [u_g, u_b]
    
    and returns:

        - Π:{Array,2} the transition matrix Π_z'ze'e
        - An entry of this matix should be read as:
            - Π_z'ze'e[z',e'] = probability of transitioning from state (z,e) to state (z',e')
    """
function trans_mat(d_z, d_u, u)

    d_z = ( d_z .- 1 )./d_z

    # transition probabilities between states: [[π_gg, π_gb][π_bg, π_bb]]
    Π_z = [d_z[1] 1-d_z[1]; 1-d_z[2] d_z[2]]

    # transition probabilities
    Π = zeros(4,4)
    d1 = Diagonal( (d_u .- 1) ./ d_u )
    Π[3:4, 3:4] = d1 + (d1 .* Diagonal([ 0.75, 1.25]))[:,end:-1:1]
    Π[1:2, 3:4] = 1 .- Π[3:4, 3:4] 
    Π[3:4, 1:2] = (u .- Π[3:4, 3:4] .* u')./(1 .- u')
    Π[1:2, 1:2] = 1 .- Π[3:4, 1:2]

    return (Π .* repeat(Π_z', 2,2) , Π_z)
end


# Set up the primitives
@with_kw mutable struct Primitives
    # Parameters of the model

    # Stochastic processes
    nZ     ::Int64          = 2                                       # number of states
    nE     ::Int64          = 2                                       # number of employment status
    d_z    ::Array{Float64} = [8, 8]                                  # Duration of states [d_g, d_b]
    d_u    ::Array{Float64} = [1.5, 2.5]                              # Duration of unemployment in each state [d_unemp_g, d_unemp_b]
    u      ::Array{Float64} = [0.04, 0.1]                             # Fraction of people unemployed in each state [u_g, u_b]
    z_val  ::Array{Float64} = [1.01, 0.99]                            # aggregate technology shocks
    e_val  ::Array{Int64}   = [1, 0]                                  # employment status

    # # # Preferences                  
    β      ::Float64          = 0.99                                   # Discount factor
    util   ::Function         = (c) -> log(c)                          # Utility function

    # # # Production
    α      ::Float64          = 0.36                                   # Capital labor ratio 
    y      ::Function         = (z, k, l) -> z * k^α * l^(1-α)         # Production function
    δ      ::Float64          = 0.025                                  # Capital depreciation rate
    ē      ::Float64          = 0.3271                                 # Labor efficiency (per hour worked)
    
    # # # Initial Conditions
    nK = 100                                                           # Number of grid points for capital
    # TODO: Derive the initial conditions from the steady state 
    # K_ss   ::Float64          = ( α  /(1 / β + δ -1 ) )^(1/(1-α))*L_ss # Steady state capital
    k_min  ::Float64          = 10.0                                   # Minimum capital
    k_max  ::Float64          = 15.0                                   # Maximum capital
    k_grid ::Array{Float64,1} = range(k_min, length = nK, stop = k_max)# Capital grid

    #L_g    ::Float64          = 1 - u[1]                               # Aggregate labor from the good state
    #L_b    ::Float64          = 1 - u[2]                               # Aggregate labor from the bad state
    #π      ::Float64          = (d_z[1] - 1) / d_z[1]                  # Long-run probability of being in good state
    #L_ss   ::Float64          = L_ss = π*L_g + (1-π)*L_b
    K_ss   ::Float64          = 11.55 #(α/((1/β) + δ - 1))^(1/(1-α))*L_ss

    K_min  ::Float64          = floor(K_ss)
    K_max  ::Float64          = 15.0
    K_grid ::Array{Float64,1} = range(K_min, length = nK, stop = K_max)# Aggregate Capital grid

    T      ::Int              = 10000                                  # Number of periods
    N      ::Int              = 5000                                   # Number of agents
    
    # First order conditions of the firm
    w      ::Function         = (K, L, z) -> (1 - α)*z*(K/L)^α
    r      ::Function         = (K, L, z) -> α*z*(K/L)^(α-1)                        # Wage rate
    
    
    # Conjecture of functional form for h₁:
    # The congecture will be a log linear function recieves the following input:
    #   - z::Int64 the technology shock
    #   - a::Array{Float64} log linear coefficients in case of good productivity shock
    #   - b::Array{Float64} log linear coefficients in case of bad productivity shock
    k_forecast ::Function     = (z, a, b, k_last) -> ( z == 1 ) ? exp(a[1]+a[2]*k_last) : exp(a[1]+a[2]*k_last)
    
end

function generate_shocks(prim)
    # TODO: This part is dependent on being two states with symmetric distributions, should be generalized
    
    @unpack N, T, d_z, d_u, u = prim
    
    # Transition Matrices
    Π, Π_z  = trans_mat(d_z, d_u, u)

    z_seq  ::Array{Int64, 1}  = vcat(1, zeros(T-1+1000))                 # Technology shocks    
    for t ∈ 2:length(z_seq)
        temp = rand(1)
        z_seq[t] = ( temp[1] < Π_z[Int(z_seq[t-1])] ) ? 1 : 2            # Generate the sequence of shocks
    end

    ℇ ::Array{Int64, 2} = zeros(N, T+1000)
    ℇ[:,1] .= 1                           # Agent's employment status
    for t ∈ 2:T
        z_last = z_seq[t-1]
        z_now  = z_seq[t]
        for n ∈ 1:N
            temp = rand(1)[1]
            e_last = ℇ[n, t-1]
            ind_1 = 2e_last + z_last
            prob_emp = Π[ind_1, z_now]
            ℇ[n,t] = ( temp < prob_emp ) ? 1 : 0
        end
    end

    # Remove the initial 1000 periods
    ℇ = ℇ[:,1001:end]
    z_seq = z_seq[1001:end]

    return (Π, Π_z,  z_seq, ℇ)
end

# Set structure of the model regarding stochastic shocks
struct Shocks
    Π      ::Array{Float64,2}                                          # Transition matrix Π_z'ze'e
    Π_z    ::Array{Float64,2}                                          # Transition matrix Π_z'z
    z_seq  ::Array{Int64, 1}                                           # Technology shocks
    ℇ      ::Array{Int64, 2}                                           # Agent's employment status
end

# Structure to hold the results
mutable struct Results
    val_fun ::Array{Float64, 2} # Value function
    pol_fun ::Array{Float64, 2} # Policy function
    a       ::Array{Float64}    # log linear coefficients in case of good productivity shock
    b       ::Array{Float64}    # log linear coefficients in case of bad productivity shock 
end

# Function to initialize the model
function Initialize()

    # Initialize the primitives
    prim = Primitives()

    # Initialize the shocks
    Π, Π_z = trans_mat(prim.d_z, prim.d_u, prim.u)                                    # Transition matrix

    # Initialize the results
    # Initialize the value function and the policy function
    val_fun = zeros(prim.nK,  prim.nE + prim.nZ)
    pol_fun = copy(val_fun)

    # Initialize the regression coefficients
    a = [0.095, 0.999] 
    b = [0.085, 0.999]

    # Return the primitives and results
    res = Results(val_fun, pol_fun, a, b)

    Π, Π_z, z_seq, ℇ = generate_shocks( prim )

    shocks = Shocks(Π, Π_z, z_seq, ℇ)

    return (prim, res, shocks)
end

# Solve consumer's problem 
function Bellman(prim::Primitives, res::Results)

    # retrieve relevant primitives and results
    @unpack k_grid, nK, nZ, nE, ē, w, r, k_forecast = prim
    @unpack a, b, pol_fun, val_fun = res



end

# Outer-most function that iterates to convergence
function SolveModel(; tol = 1e-2, err = 100, I = 1)
    
    # initialize model primitives, output, and shocks
    prim, res, shocks = Initialize()
    
    # loop until err ≤ tol 
    while err > tol

        # save last guess of a and b
        a₀ = res.a 
        b₀ = res.b

        # given current coefficients, solve consumer problem
        res = Bellman(prim, res)

        # given consumers' policy functions, simulate time series
        res, V = Simulation(prim, res, shocks)

        # Re-estimate parameters
        res = RegCapital(V)

        # update error and test convergence
        err = maximum(abs.([res.a; res.b] - [a₀; b₀]))

    end

end