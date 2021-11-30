# Load packages
using LinearAlgebra, Parameters

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
        return σ, nothing 
    end

end

# segment_data by market for demand estimation
function segment_data(model, market)

    # Get market id column
    market_id_col = model.market_id

    # Filter data by market
    data = model.inv_dem_est[model.inv_dem_est[!, market_id_col] .== market, :]
    
     # Get the observed market shares
     S = data.share
     # Get the observed prices
     P = data.price
     # Get the income levels
     Y = model.Y
     # Get the inital guess for the inverse demand
     δ = data.δ

     return S, P, Y, δ
end

# Model Structures
# Primitives
@with_kw struct Primitives
    λₚ_range :: Array{Float64} = [0, 1]
end


# Model
mutable struct Model 

    # Parameters
    parameters      ::Primitives                # Parameters of the model

    # Data
    # Unstacked data
    market_id       :: Any                      # Market id column
    product_id      :: Any                      # Product id column
    X               :: Array{Float64, 2}        # Matrix of covariates
    Z               :: Array{Float64, 2}        # Matrix of instruments
    Y               :: Array{Float64, 2}        # Matrix of simulated data
    inv_dem_est     :: DataFrame                # DataFrame of for demand estimation
    # 
end



# Demand inverter
function inverse_demand(model::Model, λₚ::Float64, market; method::String="Newton", max_iter = Inf)

    # Check the method
    valid_methods = ["Newton", "Contraction Mapping"]
    @assert (method ∈ valid_methods)

   # Get the data
    S, P, Y, δ = segment_data(model, market)

    # Compute the matrix μ[i,j] = λₚ * Y[i] * P[j]
    μ = λₚ * repeat(Y', length(S), 1) .* P

    # Compute the inverse demand
    
    # Initial guess for the inverse demand
    δ₀ = copy( δ )
    δ₁ = copy( δ )
    # Initialize the iteration counter
    iter = 0
    # Initialize the error
    err = 100
    eval_jacobian = false

    ε = 1e-12
    ε₁ = ( method == "Newton" ) ? 1e-8 : -Inf

    # Iterate until convergence

    err_list = []
    method_flag = "Contraction Mapping"
    θ = 1
    while (err > ε) && (iter < max_iter)
        # Compute the choice probability and the Jacobian
        if (method == "Newton") && (err < ε₁)
            eval_jacobian = true
            method_flag = "Newton"
        end

        σ, Δ = choice_probability(δ₀, μ, eval_jacobian=eval_jacobian)

        # Compute the inverse demand
        if (method == "Newton") && (err < ε₁)
            δ₁ = δ₀ + θ * Δ * (log.(S) - log.(σ))
        else
            δ₁ = δ₀ + log.(S) - log.(σ)
        end
        
        # Update the error
        err = maximum( abs.(δ₁ - δ₀) )
        push!(err_list, err)
        # Update the inverse demand
        δ₀ = copy(δ₁)
        
        # Update the iteration counter
        iter = iter + 1
        if iter % 1000 == 0
            println("Iteration = $iter, Method = $method_flag , error = $err, tolerance = $ε, error > tolerance = $(err > ε), θ = $θ")     
        end

        # # Update θ
        # θ = 7000 * ((12 - log10(1 / err))/2)

    end
    println("Iteration = $iter, Method = $method_flag, error = $err, tolerance = $ε, error > tolerance = $(err > ε), θ = $θ")     
    # Update the inverse demand in the model
    # for j in 1:length(J)
    #     model.δ[year][J[j]] = δ₁[j]
    # end
    println("Inverse demand = $δ₁")
    return err_list
end


