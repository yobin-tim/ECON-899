using Distributed, SharedArrays, Interpolations, Optim, Statistics

workers()
addprocs(1)

@Distributed.everywhere include("./HelpfulFunctions.jl")

Krusell_Smith()

