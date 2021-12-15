###############################################################
## This code is a julia version from inventory_problem_s2021.ox
###############################################################

using DataFrames, LinearAlgebra, CSV, Plots, Optim, Latexify

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

### Q1
## CCP mapping
function ccp(vP, lambda)
    
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

    return vEVp, aP

end

eps = 10^(-10);

err = 100

vP = fill(0.5, S) 

while err > eps 

    vP0 = vP

    tmp = ccp(vP0, lambda)

    vP = tmp[2]

    err = norm(vP0 - vP, Inf) # ox uses infinity norm

    println(err)

    if err < eps

        global vEV = tmp[1]
        global vP = tmp[2]
 
    end

end

##########################
## Q1 Output
##########################
latexify(hcat(mS, vEV))

#######################################################
## Q2
#######################################################
mF = mF0 .* (1 .- vP) + mF1 .* vP

vF = fill(1/S, S)

err = 100

while err > 10^(-10)
    
    vF0 = vF

    vF = mF' * vF0

    err = norm(vF0 - vF, Inf)

    println(err)

end

vCF = cumsum(vF)

mCF = cumsum(mF', dims = 1)'

mCF0 = cumsum(mF0', dims = 1)'

mCF1 = cumsum(mF1', dims = 1)'

## Simulation 
Sim = 5000

mSim = zeros(Sim, size(mS, 2))
    
vT = collect(Int32, 1:1:Sim)

u = rand(Float64, 1)

sid = floor(Int, sum(u.>vCF) + 1)

mSim[1,:] = mS[sid, :]

vY = zeros(Int32, Sim, 1)

vSid = zeros(Int32, Sim, 1)

vSid[1] = sid

for s = 1:Sim
    
    u = rand(Float64, 1)

    vY[s] = 1*(u .< vP[sid])[1]

    u = rand(Float64,1)

    sid = sum(u .> mCF0[sid, :], dims = 1) * (1 - vY[s]) + sum(u .> mCF1[sid, :], dims = 1) * vY[s]

    sid = floor(Int, sid[1] + 1)

    if s < Sim

        global mSim[s+1, :] = mS[sid, :]

        global vSid[s+1] = sid

    end

end
    
eta = 10^(-3)

vSid2 = vec(vSid)

vFhat = zeros(36)

vNhat = zeros(36)

vPhat = zeros(36)

for i = 1:36
    
    vFhat[i] = mean(vSid2 .== i)

    vNhat[i] = sum(vSid2 .== i)

end

for i = 1:36
    
    tmp = vNhat[i]
    
    if tmp < 1

        tmp = 1

    end
        
    vPhat[i] = (vY'* (vSid2 .== i))[1]./tmp

    if vPhat[i] < eta

        vPhat[i] = eta

    elseif vPhat[i] > 1 -eta

        vPhat[i] = 1-eta

    end

end

vEV2 = ccp(vPhat)[1]

## Q2 output
latexify(hcat(mS, vEV, vEV2))

### Likelihood function

function lfunc_ccp(lambda)
    
    vCCP = ccp(vPhat, lambda)[2]

    vCCP = vCCP[vSid]

    vL = vCCP .* vY + (1 .- vCCP) .* (1 .- vY)

    LLF = -sum(log.(vL))

end

function lfunc_nfxp(param)

    global lambda = param 
         
    global vCCP = vPhat

    global eps = 10^(-10)

    global vEV = fill(0, S) 

    global err = 100

    while err > eps 

        global vEV0 = vEV

        global vCCP0 = vCCP

        global tmp = ccp(vCCP0, lambda)

        global vEV = tmp[1]

        global err = norm(vEV - vEV0, Inf) 

        if err < eps

            global vEV = tmp[1]
            
            global vCCP = tmp[2]
 
        end

    end

    global vCCP = vCCP[vSid]

    global vL = vCCP .* vY + (1 .- vCCP) .* (1 .- vY)

    global LLF = -sum(log.(vL))

end

grid = collect(range(-10, length = 101, stop = 0))

ll = zeros(length(grid))

for i = 1:length(grid)
    
    ll[i] = lfunc_nfxp(grid[i])

end


lambda = optimize(lfunc_ccp, -10, 0).minimizer

lambda = optimize(lfunc_nfxp, -10, 0).minimizer

plot(grid, ll,
     title = "Grid serach",
     xlabel = "Lambda",
     ylabel = "Likelihood",
     color =:black,
     legend = false,
     lw = 2)

savefig("../document/figures/Q4_grid.pdf")
