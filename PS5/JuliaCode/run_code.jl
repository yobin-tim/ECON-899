using LinearAlgebra, Parameters, Plots, LaTeXStrings, Latexify
theme(:juno)
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

#include("../JuliaCode/krusell_smith.jl")
include("PS5/JuliaCode/krusell_smith.jl")

prim, res, shocks = Initialize()
V_iterate(prim, res, shocks)
V = Simulation(prim, res, shocks)

