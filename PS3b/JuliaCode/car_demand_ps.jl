#########################################################################################
## convert from ox to julia
#########################################################################################

## Load this part for each questions ###################
using Parameters, StatFiles, DataFrames, Statistics, LinearAlgebra, Plots, Optim

include("./blp_func_ps.jl")

vYear, vShare, vDelta, aProductID, mX, mZ, mEta, mIV = load_data();
## ######################################################

## Q1-1: Contraction mapping

vDelta_tmp = vDelta;

t = 1;

aJac = 0;

eps0 = 10^-12;

err_list = [];

err = 100;

iter = 0;

vParam = 0.6;

mMu = value(vParam, t);

rowid = aProductID[t];

f = 100;

while err > eps0
    
    tmp = demand(mMu, aJac, vDelta_tmp, t)

    vShat = tmp[1];

    f = log.(vShare[rowid]) - log.(vShat);

    vDelta_tmp[rowid]= vDelta_tmp[rowid] + f;

    err = norm(f)
    
    push!(err_list, err)

    iter += 1

end

x = 1:iter;

plot(x, err_list,
     title = "Contraction Mapping",
     xlabel = "Iteration",
     ylabel = "Norm",
     color =:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/Q1_contraction.pdf")

## Please shut down once because of the issue of storing data. 
## Q1-2 Combination of contraction mapping and Newton

vDelta_tmp = vDelta;

vParam = 0.6;

t = 1;

mMu = value(vParam, t)

rowid = aProductID[t]

err_list = []

err = 100

iter = 0

tmp = 0

converge = 0

f = 100

while converge == 0 

    if err > 1

        tmp = demand(mMu, 0, vDelta_tmp, 1)

        vShat = tmp[1];

        f = log.(vShare[rowid]) - log.(vShat);

        vDelta_tmp[rowid] = vDelta_tmp[rowid] + f;
    

    else
        
        tmp = demand(mMu, 1, vDelta_tmp, t)

        vShat = tmp[1];
        
        mJacobian = tmp[2];

        f = log.(vShare[rowid]) - log.(vShat);

        vDelta_tmp[rowid] = vDelta_tmp[rowid] + inv(mJacobian./vShat)*f;
        
    end

    err = norm(f)

    push!(err_list, err)

    iter += 1
    
    if err < 10^(-12)

        converge = 1

    end

end

x = 1:iter;

plot(x, err_list,
     title = "Contraction Mapping and Newton",
     xlabel = "Iteration",
     ylabel = "Norm",
     color =:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/Q1_combination.pdf")


## Please shut down once and read the head code
## Q2-Grid serch

A = inv(mIV'*mIV);

grid = collect(range(0.0, length = 11, stop = 1.0))

tmp = zeros(length(grid))

for i = 1:length(grid)
    
    tmp[i] = gmm_obj(grid[i])

end

plot(grid, tmp,
     title = "The GMM objective function",
     xlabel = "Lambda",
     ylabel = "",
     color =:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/Q2.pdf")

## Q3
## I could not use BFGS.  
lambda = optimize(gmm_obj, 0, 1).minimizer # 0.619

vDelta_init = vDelta;

res = Initialize();

inverse(res, vDelta_init, lambda)

vLParam=ivreg(res.est,mX,mIV,A);  

vXi=res.est-mX*vLParam;

mG = (vXi.*mIV);

tmp = mG .- mean(mG, dims= 1);

A=inv(tmp'*tmp);

lambda = optimize(gmm_obj, 0, 1).minimizer # 0.562
