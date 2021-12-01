#==
    This file conducts the analyses for JF's PS3
==#

# Load the necessary packages
using StatFiles, DataFrames

# Indlude the functions
include("./functions.jl")
include("./manipulate_data.jl")
include("./aux_vars.jl")

car_data, instruments, income = load_data("./PS3b/data/")
car_data, instruments, income = load_data("../data/")

# Parameters
parameters = [0.6]

# Model
model = construct_model(model_specs, car_data, instruments, income)

market, λₚ = 1985, 0.6


S, P, Y, δ = segment_data(model, market)

# Compute the matrix μ[i,j] = λₚ * Y[i] * P[j]
μ = λₚ * repeat(Y', length(S), 1) .* P

δ₀ = copy( δ )
δ₁ = copy( δ )

err, iter,err_list = 100, 0, []
while (err > 1) && (iter < 500)
    σ, Δ = choice_probability(δ₀, μ, eval_jacobian=false)

    δ₁ = δ₀ + log.(S) - log.(σ)

    # Update the error
    err = maximum( abs.(δ₁ - δ₀) )
    push!(err_list, err)
    # Update the inverse demand
    δ₀ = copy(δ₁)
    
    # Update the iteration counter
    iter = iter + 1

end

σ, Δ = choice_probability(δ₀, μ, eval_jacobian=true)

δ₁ = δ₀ +  inv(Δ) * (log.(S) - log.(σ))

# Update the error
err = maximum( abs.(δ₁ - δ₀) )
push!(err_list, err)
# Update the inverse demand
δ₀ = copy(δ₁)
