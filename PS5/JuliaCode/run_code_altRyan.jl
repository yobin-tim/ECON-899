using LinearAlgebra, Parameters, Plots, LaTeXStrings, Latexify
theme(:juno)
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

include("krusell_smith_altRyan.jl")
#include("krusell_smith.jl")
SolveModel()
