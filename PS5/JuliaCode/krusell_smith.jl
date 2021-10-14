using LinearAlgebra, Parameters, Interpolations

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
    u      ::Array{Float64} = [0.04, 0.1]                             # Fraction of people unemployed inn each state [u_g, u_b]
    z_val  ::Array{Float64} = [1.01, 0.99]                            # aggregate technology shocks
    e_val  ::Array{Int64}   = [1, 0]                                  # employment status
    e_bar  ::Float64        = 0.3271                                  # labor efficiency per unit of time worked)
    L_vals ::Array{Float64} = e_bar .* [1 - u[1] , 1 - u[2]]                              # Aggregate Labor supply
    # # # Preferences                  
    β      ::Float64          = 0.99                                   # Discount factor
    util   ::Function         = ( c ) -> log(c)                        # Utility function

    # # # Production
    α      ::Float64          = 0.36                                   # Capital labor ratio 
    y      ::Function         = ( z, k, l) -> z * k^α * l^(1-α)        # Production function
    
    # # # Initial Conditions
    # TODO: Derive the initial conditions from the steady state 
    # K_ss   ::Float64          = ( α  /(1 / β + δ -1 ) )^(1/(1-α))*L_ss # Steady state capital
    k_min  ::Float64          = 0.0                                   # Minimum capital
    k_max  ::Float64          = 20.0                                   # Maximum capital
    nk     ::Int64            = 21                                                           # Number of grid points for capital
    k_grid ::Array{Float64,1} = range(k_min, length = nk, stop = k_max)# Capital grid

    K_min  ::Float64          = 10.0
    K_max  ::Float64          = 15.0
    nK     ::Int64            = 11                                      # Number of grid points for capital
    K_grid ::Array{Float64,1} = range(K_min, length = nK, stop = K_max)# Aggregate Capital grid

    T      ::Int              = 10000                                  # Number of periods
    N      ::Int              = 5000                                   # Number of agents
    
    # First order conditions of the firm
    w_eq   ::Function         = ( K, L, z ) -> (1 - α)*z*(K/L)^α
    r_eq   ::Function         = ( K, L, z ) -> α*z*(K/L)^(α-1)                        # Wage rate
    
    
    # Conjecture of functional form for h₁:
    # The congecture will be a log linear function recieves the following input:
    #   - z::Int64 the technology shock
    #   - a::Array{Float64} log linear coefficients in case of good productivity shock
    #   - b::Array{Float64} log linear coefficients in case of bad productivity shock
    k_forecast ::Function     = (z, a, b, k_last) -> ( z == 1 ) ? exp(a[1]+a[2]*k_last) : exp(a[1]+a[2]*k_last)
    
end

function generate_shocks( prim )
    # TODO: This part is dependent on being two states with symmetric distributions, should be generalized
    
    @unpack N, T, d_z, d_u, u = prim
    
    # Transition Matrices
    Π, Π_z  = trans_mat(d_z, d_u, u)

    z_seq  ::Array{Int64, 1}  = vcat(1, zeros(T-1+1000))                    # Technology shocks    
    for t ∈ 2:length(z_seq)
        temp = rand(1)[1]
        z_seq[t] = ( temp < Π_z[ Int(z_seq[t-1])] ) ? 1 : 2          # Generate the sequence of shocks
    end

    ℇ ::Array{Int64, 2} = zeros(N, T+1000)
    ℇ[:,1] .= 1                           # Agent's employment status
    for t ∈ 2:T+1000
        z_last = z_seq[t-1]
        z_now  = z_seq[t]
        for n ∈ 1:N
            temp = rand(1)[1]
            e_last = ℇ[n, t-1]
            ind_1 = 2(1 - e_last) + z_last
            prob_emp = Π[ind_1, z_now]
            ℇ[n,t] = ( temp < prob_emp ) ? 1 : 0
        end
    end

    # Remove the initial 1000 periods
    # ℇ = ℇ[:,1001:end]
    # z_seq = z_seq[1001:end]

    return (Π, Π_z,  z_seq, ℇ)
end

# Set structure of the model regarding stochastic shocks
struct Shocks
    Π      ::Array{Float64,2}                                          # Transition matrix Π_z'ze'e
    Π_z    ::Array{Float64,2}                                          # Transition matrix Π_z'z
    z_seq  ::Array{Int64, 1}                                            # Technology shocks
    ℇ      ::Array{Int64, 2}                                           # Agent's employment status
end

# Structure to hold the results
mutable struct Resutls
    # TODO: Generalize sizes
    # We are going to define the val_fun and pol_fun as 4-dimenstional objects
    # v[:,:,z,e] gives the value functon for all posiible (k,K) combiantions for a particular (z,e) combination
    val_fun ::Array{Float64, 4} # Value function
    pol_fun ::Array{Float64, 4} # Policy function
    # we are also going to generate an iteration object for each of the functons
    # these will actually be a collection of interpolation objects for each combiantion (z,e)
    # We will store this objects in a dictionary, the idea is that the dictionary keys are (i,j) 
    # the index of z and e this is convenient for accesing purpuses later on
    val_fun_interp ::Dict
    pol_fun_interp ::Dict
    a       ::Array{Float64}    # log linear coefficients in case of good productivity shock
    b       ::Array{Float64}    # log linear coefficients in case of bad productivity shock 
end



# Function to initialize the model
function Initialize()
    # Initialize the primitives
    prim = Primitives()
    @unpack k_grid, K_grid = prim
    # Initialize the shocks
    Π, Π_z = trans_mat(prim.d_z, prim.d_u, prim.u)                                    # Transition matrix
    # Initialize the results
    # Initialize the value function and the policy function
    val_fun = zeros(prim.nk, prim.nK, prim.nZ, prim.nE)
    pol_fun = copy(val_fun)
    val_fun_interp = Dict()
    pol_fun_interp = Dict()
    # TODO: Generalize sizes
    for i ∈ 1:prim.nZ
        for j ∈ 1:prim.nE
            val_fun_interp[(i,j)] = interpolate( (k_grid, K_grid) , val_fun[:,:, i, j], Gridded(Linear() ))
            pol_fun_interp[(i,j)] = interpolate( (k_grid, K_grid) , pol_fun[:,:, i, j], Gridded(Linear() ))
        end
    end
    # Initialize the regression coefficients
    a = [0.095, 0.999] 
    b = [0.085, 0.999]
    # Return the primitives and results
    
    Π, Π_z, z_seq, ℇ = generate_shocks( prim )
    
    
    shocks = Shocks(Π, Π_z, z_seq, ℇ)
    res = Resutls(val_fun, pol_fun, val_fun_interp, pol_fun_interp, a, b)
    
    return (prim, res, shocks)
end

prim, res, shocks = Initialize()
# Solve the conumers problem 
# First we define the Bellman operator 
function Bellman( prim, res, shocks )
    
    @unpack nK, nZ, nE, k_forecast, z_vals, k_grid, K_gird, L_vals, w_eq, r_eq, e_bar, util = prim
    @unpack val_fun, pol_fun, val_fun_interp, a, b = res
    @unpack Π = shocks

    # Iterate over the stochastic states
    for i_z ∈ 1:nZ
        z = z_vals[i_z]
        L = L_vals[i_z]
        # e = e_values[i_e]
        # Iterate over all possible combinations of capital holdings and aggregate capital
        for i_K ∈ 1:nK
            K = prim.K_grid[i_K]
            k = prim.k_grid[i_k]
            # Consumers forecast next period aggregate capital
            K_Next = k_forecast(z, a, b, K )
            # Itertate over employment states
            w = w_eq(K, L, z)
            r = r_eq(K, L, z)
            for i_e ∈ 1:nE   
                # Next we need to get wich row of the markov matrix are we
                row = 2*i_z + i_e
                ##### HERE
                # Iterate over all posible values for capital 
                for i_k ∈ 1:nK
                    # Calculte budget given current capital
                    budget = r * k + w * e_bar * e + (1 - δ)*k
                    # Iterate over all posible next period value for capital
                    for i_k_next ∈ 1:nK
                        k_next = k_grid[i_k_next]
                        # Calculate consumption
                        c = budget - k_next
                        # Calculate the Utility from consumption
                        utility = util(c)
                        # Calculate the continuation value
                        # For this we will use the interpolated version of the value function
                        # since K_tomorrow may not be in the grid
                        # To calculate the expeted vaue we will multiply the row of the trnasition matrix
                        # to a column vector containing the values of each of the posible values of the state of the world
                        # given next period aggregate capital and individual capital choice
                        exp_val_next = Π[row, :]
                    end
                end
            end
        end
    end

end # Bellman 

# Next we will do value function iteration to solve the conusmer problem
function Solve_Bellman( prim, res )

end

using Interpolations

x = 0:0.01:5
y = 0:0.01:5

z = zeros(length(x), length(y))

for i ∈ 1:length(x), j ∈ 1:length(y)
    z[i,j] = x[i]^2 + y[j]^2
end

x_data = 0:0.1:5
y_data = 0:0.1:5


z_data = zeros(length(x_data), length(y_data))
for i ∈ 1:length(x_data), j ∈ 1:length(y_data)
    z_data[i,j] = x_data[i]^2 + y_data[j]^2
end


interp = interpolate((x_data, y_data), z_data, Gridded(Linear()))



interp.( 4.99999, 4.99999 )

interp_z = zeros(length(x), length(y))
for i ∈ 1:length(x), j ∈ 1:length(y)
    interp_z[i,j] = interp(x[i], y[j])
end
# Sup norm interpdata - actual function 
maximum( abs.( interp_z - z  ) )
