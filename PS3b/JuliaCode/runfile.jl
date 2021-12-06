#==
    This file conducts the analyses for JF's PS3
==#

# Load the necessary packages
using StatFiles, DataFrames, Plots, LaTeXStrings

# Indlude the functions
include("./functions.jl")
include("./manipulate_data.jl")
include("./aux_vars.jl")

car_data, instruments, income = load_data("./PS3b/data/")
#car_data, instruments, income = load_data("../data/")

# Parameters
parameters = [0.6]

# Model
model = construct_model(model_specs, car_data, instruments, income)

market, λₚ = 1985, parameters[1]

# Testing the GMM part
    l = 0:0.01:0.5
    data = gmm.(Ref(model), l)

    plot(l, data)

    model = construct_model(model_specs, car_data, instruments, income)

    inverse_demand(model, λₚ, market; method = "Contraction Mapping", max_iter = 100)

    S, P, Y, δ = segment_data(model, market)

    μ = λₚ * repeat(Y', length(S), 1) .* P

    δ₀ = copy( δ )
    δ₁ = copy( δ )

    σ, Δ = choice_probability(δ₀, μ, eval_jacobian=true)

    δ₁ = δ₀ + inv(Δ) * (log.(σ) - log.(S))

    err = maximum( abs.(δ₁ - δ₀) )

    δ₀ = copy(δ₁)

###############################################################################
####                            Problem 1
###############################################################################
model = construct_model(model_specs, car_data, instruments, income)
markets = unique(car_data.Year)

err_list = Dict()
for market in markets
    err_list[market] = inverse_demand(model, λₚ, market)
end
err_list

p = plot([err_list[market][2:end] for market in markets], label="")
scatter(err_list[2000][2:end])



#inverse_demand(model, λₚ, 1992)

###############################################################################
####                            Problem 2
###############################################################################
l = 0:0.01:1
data = gmm.(Ref(model), l)
    plot(l, data,title="GMM Objective Function", xlabel=L"\lambda_{p}", legend=false)
    savefig("./PS3b/Document/Figures/Problem2.png")


###############################################################################
####                            Problem 3
###############################################################################
    λhat_GMM=TwoStage_gmm(model)
