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

vP = fill(0.5, S) 
