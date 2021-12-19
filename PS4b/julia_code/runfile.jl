#==
    This file conducts the analyses for JF's PS4
==#

#==
        Housekeeping
    ---------------------
==#


# Load the necessary packages
using DataFrames, LinearAlgebra, CSV, Plots, Optim, Latexify, PGFPlotsX
theme(:juno)
# Un-comment the following line when compiling the document
# theme(:vibrant)
default(fontfamily = "Computer Modern", framestyle = :box) # LaTex-style


# Indlude the functions
include("./functions.jl")

# Initialize the model primitives
prim = Primitives(λ = 0.0)

# load the simulated data
sim = (DataFrame(CSV.File("./PS4b/ox_code/PS4_simdata.csv")) 
        |> Matrix)[:, 2:end]



#==
        Problem 1
    -------------------
==#

# initial guess of P: 1/2 for all entries
P₀ = .5*ones(size(prim.S, 1), 1)

# Iterate until CCP convergence 
EV, P = CCP(prim, P₀; N = 1);

# print the fixed-point EV
fname = "./PS4b/document/Problem1.tex";
open(fname, "w") do io
    write(io, latexify(round.(EV, digits = 2)))
end;



#==
        Problem 2
    -------------------
==#

# simply calculate Phat by calculating the share of times a=1 for 
# each state 
Phat = [sum(sim[sim[:, 2] .== i, 1] .== 1) / sum(sim[:, 2] .== i) 
            for i ∈ 0:sort(unique(sim[:, 2]))[end]]

# winsorize at .001 and .999
Phat[Phat .< 0.001] .= 0.001; Phat[Phat .> 0.999] .= 0.999;

# Calculate the implied calue function 
EVhat = ExpVal(prim, Phat, prim.λ₀)[1]

# print the resulting EV, alongside the result from problem 1
fname = "./PS4b/document/Problem2.tex";
open(fname, "w") do io
    write(io, latexify(round.([EV EVhat], digits = 2)))
end;



#==
        Problem 4
    -------------------
==#

# plot log-likelihood function from -10 to 0
grid   = range(-10, 0, length = 1000)
LLfunc = [-LL(sim[:, 1], sim[:, 2], grid[i], Phat) for i ∈ 1:length(grid)] 

plot(grid, LLfunc,
    xlabel = "λ",
    title = "Log-likelihood Function (negative)",
    legend = nothing
)
savefig("./PS4b/document/figures/Problem4.pdf")

# find the λ that maximizes LL
λ = optimize(x -> -LL(sim[:, 1], sim[:, 2], x, Phat), -10.0, 0.0).minimizer

# print the result
fname = "./PS4b/document/Problem4.tex";
open(fname, "w") do io
    str = @sprintf "\$\\hat{\\lambda}=%1.3f\$" λ
    write(io, str)
end;