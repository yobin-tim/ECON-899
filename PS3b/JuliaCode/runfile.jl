#==
    This file conducts the analyses for JF's PS3
==#

# Load the necessary packages
using StatFiles, DataFrames

# Indlude the functions
include("./functions.jl")
include("./manipulate_data.jl")

car_data, instruments, income = load_data("./PS3b/data/")

# Parameters
parameters = [0.6]

# Model
model = construct_model(car_data, instruments, income, parameters)

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
r,R = size(μ) 

δ = inverse_demand_cm(δ, S, μ, 1)
δ = inverse_demand_nm(δ, S, μ, 0.000001)



σ, Δ  = choice_probability(δ, μ; eval_jacobian=true)

Δ
1/R * ((I(r) .* (σ * (1 .- σ)'))./ σ- ((1 .- I(r)) .* (σ * σ'))./ σ)

inv(Δ)

δ₁ = δ - inv(Δ) * (log.(S) - log.(σ))

err = maximum( abs.(δ₁ - δ) ) 

δ = copy(δ₁)

for i in 1:length(δ)
    println(abs.(δ₁[i] - δ[i]))
end 
# Not working yet
# inverse_demand(model, 1985; method = "Newton")

inverse_demand(model, 2005, method = "Newton")

model.δ[1985]

