# Load packages
using LinearAlgebra

# Auxiliary functions
# Choice probability function
function choice_probability(δ::Array{Float64}, μ::Array{Float64}; eval_jacobian::Bool = false)
    
    # number of individuals and choicesm
    r, R = size(μ);

    # Compute choice probabilities
    Λ = exp.( δ .+ μ  )
    Λ = Λ./(1 .+ sum(Λ, dims=1))
    σ = sum(Λ, dims=2)/R

    if eval_jacobian
        # Compute Jacobian
        Δ = 1/R * ((I(r) .* (σ * (1 .- σ)')) - ((1 .- I(r)) .* (σ * σ')))
        return σ, Δ
    else
        return σ   
    end

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
function inverse_demand_cm(δ, s, μ, tol; max_iter::Int64 = 300)

    # Initial guess for the inverse demand
    δ₀ = copy( δ )
    δ₁ = copy( δ )
    
    # Initialize the iteration counter
    iter = 0
    # Initialize the error
    err = 100

    # Iterate until convergence
    while (err > tol) && (iter < max_iter)
        # Compute the choice probability
        σ = choice_probability(δ₀, μ)

        # Compute the inverse demand
        δ₁ = δ₀ + log.(s) - log.(σ)
        
        # Update the error
        err = maximum( abs.(δ₁ - δ₀) )

        # Update the inverse demand
        δ₀ = copy(δ₁)
        
        # Update the iteration counter
        iter = iter + 1
        println("Iteration = $iter , error = $err, tolerance = $tol, error > tolerance = $(err > tol)")
    end

    return δ₁
end

# Demand inversion Newtons Method
function inverse_demand_nm(δ, s, μ, tol; max_iter::Int64 = 300)

    # Initial guess for the inverse demand
    δ₀ = copy( δ )
    δ₁ = copy( δ )
    # Initialize the iteration counter
    iter = 0
    # Initialize the error
    err = 100

    # Iterate until convergence
    while (err > tol) && (iter < max_iter)
        # Compute the choice probability and the Jacobian
        σ, Δ = choice_probability(δ₀, μ, eval_jacobian=true)

        # Compute the inverse demand
        δ₁ = δ₀ + Δ * (log.(s) - log.(σ))
        
        # Update the error
        err = maximum( abs.(δ₁ - δ₀) )

        # Update the inverse demand
        δ₀ = copy(δ₁)
        
        # Update the iteration counter
        iter = iter + 1
        println("Iteration = $iter , error = $err, tolerance = $tol, error > tolerance = $(err > tol)")
    end
    return δ₁
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

    # Compute the matrix μ[i,j] = λₚ * Y[i] * P[j]
    μ = repeat(Y', length(J), 1) .* P

    # Compute the inverse demand

    if method=="Newton"
        ε = 1e-12
        ε₁ = 1
        println("Aprox inverse demand usign Contraction Mapping")
        δ = inverse_demand_cm(δ, S, μ, ε₁)
        println("Solving inverse demand usign Newton's Method")
        δ = inverse_demand_nm(δ, S, μ, ε)
    elseif method=="CM"
        ε = 1e-12
        println("Solving inverse demand usign Contraction Mapping")
        δ = inverse_demand_cm(δ, S, μ, ε)
    else
        error("Method not implemented, implemented methods are \"Newton\" and \"CM\" ")
    end

    # Update the inverse demand in the model
    for j in 1:length(J)
        model.δ[year][J[j]] = δ[j]
    end

end


