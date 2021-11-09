using LinearAlgebra, Parameters, Plots, LaTeXStrings, Latexify
theme(:juno)
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

include("krusell_smith.jl")

SolveModel()
