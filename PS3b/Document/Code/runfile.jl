#==
    This file conducts the analyses for JF's PS3
==#

# Load the necessary packages
using StatFiles, DataFrames, Plots, LaTeXStrings, Printf
theme(:juno)
# Un-comment the following line when compiling the document
# theme(:vibrant)
default(fontfamily = "Computer Modern", framestyle = :box) # LaTex-style


# Indlude the functions
include("./functions.jl")
include("./manipulate_data.jl")
include("./aux_vars.jl")

#car_data, instruments, income = load_data("./PS3b/data/")
# car_data, instruments, income = load_data("../data/")
car_data, instruments, income = load_data("data")



###############################################################################
####                            Problem 1
###############################################################################
model = construct_model(model_specs, car_data, instruments, income)
market, λₚ = 1985, 0.6

err_list_cm = inverse_demand(model, λₚ, market; method = "Contraction Mapping", max_iter = Inf)
# No longer need to reset the model. Funcion will recalculate δ from scratch if not told otherwise
err_list_n = inverse_demand(model, λₚ, market; method = "Newton", max_iter = Inf)

plot(err_list_cm[2:end], xlabel = "Iteration", size = (800, 600),
    ylabel = "Error", title = "Error vs. Iteration", label = "Contraction Mapping")
plot!(err_list_n[2:end], line = (:dash), label = "Newton")
    #savefig("./PS3b/Document/Figures/Problem1.pdf")
    savefig("Document/Figures/Problem1.pdf")


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


grid = 0:0.1:1
#? Experimemt
# I'll try different values of the parameter and check if using previous demand estimates improves speed of convergence

model = construct_model(model_specs, car_data, instruments, income)

num_iter_recalc = []
for λ ∈ grid
    err = inverse_demand(model, λ, market; recalc_δ = true)
    push!(num_iter_recalc, length(err))
end

num_iter_no_recalc = []
for λ ∈ grid
    err = inverse_demand(model, λ, market; recalc_δ = false)
    push!(num_iter_no_recalc, length(err))
end


sum(num_iter_recalc) / length(num_iter_recalc)
sum(num_iter_no_recalc) / length(num_iter_no_recalc)

#? Experimemt
# What is faster Newton or Contraction Mapping?
using BenchmarkTools
@btime begin
    model = construct_model(model_specs, car_data, instruments, income)
    for λ ∈ grid
        err = inverse_demand(model, λ, market; recalc_δ = false, method = "Newton")
    end
end

@btime begin
    model = construct_model(model_specs, car_data, instruments, income)
    for λ ∈ grid
        err = inverse_demand(model, λ, market; recalc_δ = false, method = "Contraction Mapping")
    end
end

function ReturnData(model, grid)
    data = []
    for λ in grid
        println(λ)
        data = push!(data, gmm(model, λ))
    end
    return data
end

data = ReturnData(model, grid)

plot(grid, data, title = "GMM Objective Function", xlabel = L"\lambda_{p}", legend = false, markershape = :auto)
xticks!(grid)
#savefig("./PS3b/Document/Figures/Problem2.pdf")
savefig("Document/Figures/Problem2.pdf")

###############################################################################
####                            Problem 3
###############################################################################

λhat_GMM = TwoStage_gmm(model)

# print result
#fname = "./PS3b/Document/Problem3.tex";
fname = "Document/Problem3.tex";
open(fname, "w") do io
    str = @sprintf "\$\\lambda_p=%1.3f\$" λhat_GMM[1]
    write(io, str)
end;



#==
##Testing
opt = optimize(λ -> gmm(model, λ[1]), [.6], method = BFGS(), f_tol = 1e-5, g_tol = 1e-5)
    λ_hat = opt.minimizer[]
    ξ_hat=gmm(model, λ_hat, Return_ρ=true)
    OptimalW=inv( (model.Z .* ξ_hat)*transpose(model.Z .* ξ_hat) )
    λ_hat_2=gmm(model, λ_hat,SpecifyW=true,SpecifiedW=OptimalW)


Z = model.Z

temp = (Z' * ξ_hat) * (Z' * ξ_hat)'



det(temp .+ rand(0.9:0.001:1.1, size(temp)))

size(temp .+ rand(0.9:0.001:1.1, size(temp)))



OptimalW=inv( inv(temp .+ rand(0.9:0.001:1.1, size(temp))) )

Z

Zopt_scondstage = optimize(λ -> gmm(model, λ[1],SpecifyW=true,SpecifiedW=OptimalW),
                [λ_hat], method = BFGS(), f_tol = 1e-5, g_tol = 1e-5).minimizer

TwoStage_gmm(model)
==#
