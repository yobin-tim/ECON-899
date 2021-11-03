# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 09/07/2021
# neogrowth_script_parallel.jl
# ------------------------------------------------------------------------------
using Distributed
addprocs(6)

@everywhere using Parameters, Plots, SharedArrays
dir = "/Users/philipcoyle/Documents/School/University_of_Wisconsin/ThirdYear/Fall_2021/TA - Computation/ProblemSets/PS1/Julia"
cd(dir)
@everywhere include("neogrowth_model_parallel.jl")

## Main Code
@elapsed G, PFs = solve_model()
# plot_pfs(dir, G, PFs)
