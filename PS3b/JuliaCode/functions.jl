# Load packages
using Parameters

# Auxiliary functions
# Choice probability function
function choice_probability(δ::Array{Float64}, μ::Array{Float64}; eval_jacobian::Bool = false)
    
    R = size(μ)[2];
    Δ = exp.( δ .+ μ  )
    Δ = Δ./(1 .+ sum(Δ, dims=1))

    # if eval_jacobian
        #return Δ, jacobian(Δ, δ)
    # else
        return sum(Δ, dims=2)/R
    # end
end

# Model Structures

# Model
mutable struct Model 

    # Parameters
    λₚ      :: Float64 # Somehting

    # Data
    Years   :: Array{Int64}  # Array of years
    J       :: Dict          # Array of product indexes
    S       :: Dict          # Dictionary of observed market shares
    P       :: Dict          # Dictionary of observed prices
    Y       :: Array{Float64} # Array of simulated income levels

    # Functions
    # choice_probability :: Function(Float64, Float64) -> Float64

    # Demand
    δ       :: Dict # Dictionary of (estimated) inverse demand

end


# Demand inversion Contraction Mapping (Barry et al., 1995)
function inverse_demand_cm(δ, s, μ, tol; max_iter::Int64 = 1000)

    # Initial guess for the inverse demand
    δ₀ = copy(δ )

    # Initialize the iteration counter
    iter = 0
    # Initialize the error
    err = 100

    # Iterate until convergence
    while (err > tol) && (iter < max_iter)
        # Compute the choice probability
        σ = choice_probability(δ₀, μ)

        # Compute the inverse demand
        δ = δ₀ + log.(s) - log.(σ)
        
        # Update the error
        err = maximum( abs.(δ - δ₀) )
        
        # Update the inverse demand
        δ₀ = copy(δ)
        
        # Update the iteration counter
        iter = iter + 1
        @show iter
        @show err
    end

    return δ
end

# Demand inverter
function inverse_demand(model::Model, year::Int64; method::String="Newton")

    # Get the product indexes
    J = model.J[year]
    # Get the observed market shares
    S = [model.S[year][j] for j in J]
    # Get the observed prices
    P = [model.P[year][j] for j in J]
    # Get the income levels
    Y = model.Y
    # Get the inital guess for the inverse demand
    δ = [model.δ[year][j] for j in J]

    μ = repeat(Y', length(J), 1) .* P




    if method=="Newton"

    elseif method=="CM"

    else
        error("Method not implemented, implemented methods are \"Newton\" and \"CM\" ")
    end

end


