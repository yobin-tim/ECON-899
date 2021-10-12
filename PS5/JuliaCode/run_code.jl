using LinearAlgebra, Parameters, Plots, LaTeXStrings, Latexify
theme(:juno)
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

include("../JuliaCode/krusell_smith.jl")

prim, res = Initialize()


