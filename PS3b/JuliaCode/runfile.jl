#==
    This file conducts the analyses for JF's PS3
==#

# Load the necessary packages
using StatFiles, DataFrames, Plots, LaTeXStrings
theme(:juno)
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style


# Indlude the functions
include("./functions.jl")
include("./manipulate_data.jl")
include("./aux_vars.jl")

car_data, instruments, income = load_data("./PS3b/data/")
# car_data, instruments, income = load_data("../data/")



###############################################################################
####                            Problem 1
###############################################################################
model = construct_model(model_specs, car_data, instruments, income)
market, λₚ = 1985, 0.6

err_list_cm = inverse_demand(model, λₚ, market; method = "Contraction Mapping", max_iter = Inf)
    #Reset the model
    model = construct_model(model_specs, car_data, instruments, income)
    market, λₚ = 1985, 0.6
err_list_n = inverse_demand(model, λₚ, market; method = "Newton", max_iter = Inf)

plot(err_list_cm[50:end],50:length(err_list_cm), xlabel = "Iteration",
    ylabel = "Error", color=[:white],title = "Error vs. Iteration", label="Contraction Mapping")
plot!(err_list_n[50:end],50:length(err_list_n), color=[:blue],line=(:dash), label="Newton")
savefig("./PS3b/Document/Figures/Problem1.png")

# markets = unique(car_data.Year)

# err_list = Dict()
# for market in markets
#     err_list[market] = inverse_demand(model, λₚ, market)
# end
# err_list

# p = plot([err_list[market][2:end] for market in markets], label="")
# scatter(err_list[2000][2:end])

#inverse_demand(model, λₚ, 1992)

###############################################################################
####                            Problem 2
###############################################################################

function ReturnData(model)
    l = 0:0.01:.1
    data = []
    for λ in l
        println(λ)
        data = push!(data, gmm(model, λ))
    end
    return data
end
data=ReturnData(model)
    plot(l, data, title="GMM Objective Function", xlabel=L"\lambda_{p}", legend=false)
        savefig("./PS3b/Document/Figures/Problem2.png")

###############################################################################
####                            Problem 3
###############################################################################
 λhat_GMM=TwoStage_gmm(model)


TwoStage_gmm(model)
