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

