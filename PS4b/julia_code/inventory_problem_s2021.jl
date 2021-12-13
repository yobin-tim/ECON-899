###############################################################
## This code is a julia version from inventory_problem_s2021.ox
###############################################################

using DataFrames, LinearAlgebra, CSV

alpha = 2

lambda = -4

beta = 0.99

ibar = 8

ps = 1

pr = 4

c = 1

iI = 1

iC = 2

iP = 3

gamma = 1/2

mPi = [0.75 0.25; 0.9 0.1]

mGamma = [gamma 1-gamma; gamma 1-gamma];

vIgrid = collect(0:c:ibar)

vPgrid = [pr;ps]

vCgrid = [0;c]

S = size(vIgrid)[1] * size(vCgrid)[1] * size(vPgrid)[1] # 36

mS = DataFrame(CSV.File("../ox_code/PS4_state_space.csv")) |> Matrix

mS = mS[:,3:5]

mF0 = DataFrame(CSV.File("../ox_code/PS4_transition_a0.csv")) |> Matrix

mF0 = mF0[:,3:end]

mF1 = DataFrame(CSV.File("../ox_code/PS4_transition_a1.csv")) |> Matrix

mF1 = mF1[:,3:end]

## CCP mapping

# function value(vEV)

#     tmp = vU0 + beta * mF0 * vEV

#     tmp2 = vU1 + beta * mF1 * vEV

#     mV = hcat(tmp, tmp2)

#     return mV

# end

function ccp(vP)
    
    vU1 = alpha.*mS[:,iC] - mS[:, iP];
    
    vU0 = alpha.*mS[:,iC].*(mS[:,iI].>0)+lambda.*(mS[:,iC].>0).*(mS[:,iI].==0);

    vE1 = 0.577216 .- log.(vP);

    vE0 = 0.577216 .- log.(1 .- vP);

    mF = mF0 .* (1 .- vP) .+ mF1 .* vP;

    vEU = (1 .- vP) .* (vU0 .+ vE0) .+ vP .* (vU1 .+ vE1);

    vEVp = inv(I(size(vP,1)) .- beta .* mF)*vEU;

    vU1 = alpha.*mS[:,iC] - mS[:, iP];
    
    vU0 = alpha.*mS[:,iC].*(mS[:,iI].>0)+lambda.*(mS[:,iC].>0).*(mS[:,iI].==0);

    tmp = vU0 + beta * mF0 * vEVp

    tmp2 = vU1 + beta * mF1 * vEVp

    mV = hcat(tmp, tmp2)

    aP = (exp.(mV[:,2])) ./ (sum(exp.(mV), dims = 2));

    return mV, aP

end

eps = 10^(-10);

err = 100

vP = fill(0.5, S) 

while err > eps 

    vP0 = vP

    tmp = ccp(vP0)

    vP = tmp[2]

    err = norm(vP0 - vP, Inf) # ox uses infinity norm

    println(err)

end
