#########################################################################################
## convert from ox to julia
#########################################################################################

using Parameters, StatFiles, DataFrames, Statistics, LinearAlgebra, Plots

include("./blp_func_ps.jl")

vYear, vShare, vDelta, aProductID, mX, mZ, mEta, mIV = load_data();

A = inv(mIV' * mIV)

vParam = 0.6;

grid = collect(range(0.0, length = 11, stop = 1.0))

tmp = zeros(length(grid))

for i = 1:length(grid)
    
    tmp[i] = -gmm_obj(grid[i])[1]

end

plot(grid, -tmp,
     title = "The GMM objective function",
     xlabel = "Lambda",
     ylabel = "",
     color =:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/Q2.pdf")


## Q1-2 Combination of contraction mapping and Newton

inverse(res, vDelta_init, vParam) 
    
x = 1:iter;

plot(x, err_list,
     title = "Contraction Mapping and Newton",
     xlabel = "Iteration",
     ylabel = "Norm",
     color =:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/Q1_combination.pdf")

## Q2-Grid serch

## Q1-1: Contraction mapping

# t = 1

# vDelta = vDelta_iia

# aJac = 0

# eps0 = 1e-12

# err_list = []

# err = 100

# iter = 0

# while err > eps0
    
#     tmp = demand(mMu, aJac, vDelta, t)

#     vShat = tmp[1];

#     f = log.(vShare[rowid]) - log.(vShat);

#     vDelta[rowid]= vDelta[rowid]+f;

#     err = norm(f)
    
#     push!(err_list, err)

#     iter += 1

# end

# x = 1:iter;

# plot(x, err_list,
#      title = "Contraction Mapping",
#      xlabel = "Iteration",
#      ylabel = "Norm",
#      color =:black,
#      legend = false,
#      lw = 2)

# savefig("../Document/Figures/Q1_contraction.pdf")

