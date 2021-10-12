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
function trans_mat(d_z, d_unemp, u)
    d_z = ( d_z .- 1 )./d_z
    # transition probabilities between states: [[π_gg, π_gb][π_bg, π_bb]]
    Π_z = [d_z[1] 1-d_z[1]; 1-d_z[2] d_z[2]]
    # transition probabilities
    Π = zeros(4,4)
    d1 = Diagonal( (d_unemp .- 1) ./ d_unemp )
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
    d_z    ::Array{Float64} = [8, 8]                                  # Duration of states [d_g, d_b]
    d_u    ::Array{Float64} = [1.5, 2.5]                              # Duration of unemployment in each state [d_unemp_g, d_unemp_b]
    u      ::Array{Float64} = [0.04, 0.1]                             # Fraction of people unemployed inn each state [u_g, u_b]
    z_val  ::Array{Float64} = [1.01, 0.99]                            # aggregate technology shocks
    e_val  ::Array{Int64}   = [1, 0]                                   # aggregate productivity shocks

    # # # Preferences                  
    β      ::Float64          = 0.99                                   # Discount factor
    util   ::Function         = ( c ) -> log(c)                        # Utility function

    # # # Production
    α      ::Float64          = 0.36                                   # Capital labor ratio 
    y      ::Function         = ( z, k, l) -> z * k^α * l^(1-α)        # Production function
    
    # # # Initial Conditions
    nK = 100                                                           # Number of grid points for capital
    # TODO: Derive the initial conditions from the steady state 
    # K_ss   ::Float64          = ( α  /(1 / β + δ -1 ) )^(1/(1-α))*L_ss # Steady state capital
    k_min  ::Float64          = 10.0                                   # Minimum capital
    k_max  ::Float64          = 15.0                                   # Maximum capital
    k_grid ::Array{Float64,1} = range(k_min, length = nK, stop = k_max)# Capital grid

    K_min  ::Float64          = 11.0
    K_max  ::Float64          = 15.0
    K_grid ::Array{Float64,1} = range(K_min, length = nK, stop = K_max)# Aggregate Capital grid

    # # TODO: This part is dependent on being two states with symmetric distributions, should be generalized
    T      ::Int              = 11000                                  # Number of periods
    N      ::Int              = 5000                                   # Number of agents

    

    # Conjecture of functional form for h₁:
    # The congecture will be a log linear function recieves the following input:
    #   - z::Int64 the technology shock
    #   - a::Tuple{Float64, Float64} log linear coefficients in case of good productivity shock
    #   - b::Tuple{Float64, Float64} log linear coefficients in case of bad productivity shock
    k_forecast ::Function     = (z, a, b, k_last) -> ( z == 1 ) ? exp(a[1]+a[2]*k_last) : exp(a[1]+a[2]*k_last)
    
end

function generate_shocks( N, T, d_z, d_u, u   )
    
    @unpack  = prim
    
    z_seq  ::Array{Int64, 1}  = vcat(1, zeros(T-1))                    # Technology shocks    
    for t ∈ 2:length(shoks)
        temp = rand(1)
        z_seq[t] = ( temp < Π_z[ Int(z_seq[t-1])] ) ? 1 : 2          # Generate the sequence of shocks
    end

    ℇ = zeros(N, T)                             # Agent's employment status
    for t ∈ 1:T
        z_last = z_seq[t-1]
        z_now  = z_seq[t]
        for n ∈ 1:N
            temp = rand(1)
            e_last = ℇ[n, t-1]
            ind_1 = 2e_last + z_last
            prob_emp = Π[ind_1, z_now]
            ℇ[n,t] = ( temp < prob_emp ) ? 1 : 0
        end
    end
    return (z_seq, ℇ)
end

# Set structure of the model regarding stochastic shocks
struct shocks()
    Π      ::Matrix{Float64,2}                                          # Transition matrix Π_z'ze'e
    Π_z    ::Matrix{Float64,2}                                          # Transition matrix Π_z'z
    Π, Π_z = 
    
end

# Structure to hold the results
mutable struct Resutls
    val_fun ::Array{Float64, 2} # Value function
    pol_fun ::Array{Float64, 2} # Policy function
end

# Function to initialize the model
function Initialize()
    # Initialize the primitives
    prim = Primitives()
    # Initialize the shocks
    Π, Π_z = trans_mat(d_z, d_u, u)                                    # Transition matrix
    # Initialize the results
    # Initialize the value function and the policy function
    val_fun = zeros(prim.nK, maximum(size(prim.Π)))
    pol_fun = copy(val_fun)
    # Return the primitives and results
    res = Resutls(val_fun, pol_fun)
    return (prim, res)
end

prim, res = Initialize()