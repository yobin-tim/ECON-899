using Parameters, DelimitedFiles

# Define the primitives of the model
@with_kw mutable struct  Primitives
    N    ::Int64 = 66         # Lifespan of the agents
    J_R  ::Int64  = 46        # Retirement age
    n    ::Float64  = 0.011   # Population growth rate
    a_1  ::Float64  = 0       # Initial assets holding for newborns
    θ    ::Float64  = 0.11    # Labor income tax rate
    γ    ::Float64 = 0.42     # Utillity weight on consumption
    σ    ::Float64 = 2.0      # coefcient of relative risk aversion
    α    ::Float64  = 0.36    # Capital share in production
    δ    ::Float64  = 0.06    # Capital depreciation rate
    
    # Parameters regarding stochastic processes
    z_H  ::Float64  = 1.0     # Idiosyncratic productivity High
    z_L  ::Float64  = 0.5     # Idiosyncratic productivity Low
    p_H  ::Float64  = 0.2037  # Probability of z_H at birth
    p_L  ::Float64  = 0.7963  # Probability of z_L at birth
    # Markov transition matrix for z
    Π    ::Array{Float64,2} = [0.9261 1-0.9261;  0.9811 1- 0.9811] 

    # Functions
    utility_W::Function = (c, l) -> (c^γ * (1 - l)^γ)^(1-σ)/(1-σ) # Utility of a worker
    utility_R::Function = (c) -> c^(γ*(1-σ))/(1-σ)

    # Grids
    η::Array{Float64,1} = readdlm("Shared Repo/PS3/Data/ef.txt") # Age efficiency profile
    
end

workspace()