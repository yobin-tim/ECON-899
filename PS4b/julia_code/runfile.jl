#==
    This file conducts the analyses for JF's PS4
==#

#==
        Housekeeping
    ---------------------
==#


# Load the necessary packages
using DataFrames, LinearAlgebra, CSV, Plots, Optim, Latexify
theme(:juno)
# Un-comment the following line when compiling the document
# theme(:vibrant)
default(fontfamily = "Computer Modern", framestyle = :box) # LaTex-style


# Indlude the functions
include("./functions.jl")

# Initialize the model primitives
prim = Primitives()

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
#==
#TODO: figure out how to print vector to string
# print the fixed-point EV
fname = "./PS4b/document/Problem1.tex";
open(fname, "w") do io
    str = @sprintf "%.2f" EV
    write(io, str)
end;
==#
