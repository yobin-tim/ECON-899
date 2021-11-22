#==
    This file conducts the analyses for JF's PS2
==#

using StatFiles, DataFrames, Optim, BenchmarkTools, Latexify, CSV
# We can use BenchmarkTools for better precision. Just need to replace time
# with btime. The runtime of the overall code get longer as btime runs the
# code multiple times to reduce noise

#include("./PS2b/JuliaCode/functions.jl")
include("./functions.jl")

## load the mortgage data and sparse grid weights as a DataFrames (and
## convert weights to matrices
df = DataFrame(StatFiles.load("PS2b/data/Mortgage_performance_data.dta"))
w1 = DataFrame(CSV.File("PS2b/data/KPU_d1_l20.csv")) |> Matrix
w2 = DataFrame(CSV.File("PS2b/data/KPU_d2_l20.csv")) |> Matrix

# Use this if you are loading data from the root folder.
df = DataFrame(StatFiles.load("../data/Mortgage_performance_data.dta"))
w1 = DataFrame(CSV.File("../data/KPU_d1_l20.csv")) |> Matrix
w2 = DataFrame(CSV.File("../data/KPU_d2_l20.csv")) |> Matrix

#df = DataFrame(CSV.File("C:/Users/ryana/OneDrive/Documents/School/PhD Economics/Research/GitHub/ECON-899/PS1b/data/mortgage.csv"))

## Separate data into independent variable matrices X and Z and
## dependent variable vector Y
X = df[!, [:score_0, :rate_spread, :i_large_loan, :i_medium_loan,
    :i_refinance, :age_r, :cltv, :dti, :cu,
    :first_mort_r, :i_FHA, :i_open_year2,
    :i_open_year3, :i_open_year4, :i_open_year5]] |> Matrix;

Z = df[!, [:score_0, :score_1, :score_2]] |> Matrix;

Y = df[!, :duration]; #|> Matrix


for name in names(df)
    println(name)
end
## 1. Evaluate log-likelihood using the quadrature method
println("See QuadLL() function")

## 2. Evaluate simulated log-likelihood function using GHK
println("see GHKLL() function")

## 3. Evaluate simulated log-likelihood function using accept/reject
println("see AcceptRejectLL() function")


## 4. Compare predicted choice probabilities for each of the above methods
θ₀ = [0, -1, -1, 0 * ones(size(X, 2), 1), 0.3 * ones(size(Z, 2), 1), 0.5]

ll_quad=QuadLL2(Y, X, Z, w1, w2, θ₀)
ll_ghk=GHKLL(Y, X, Z, θ₀)
ll_ar=AcceptRejectLL(Y, X, Z, θ₀)

## 5. Maximize quadrature log-likelihood function using BFGS
# TODO: This method gives domain issues with the logs in QuadLL
θ₀ = vcat([0, -1, -1], 0 * ones(size(X, 2), 1), 0.3 * ones(size(Z, 2), 1), [0.5])

θ = optimize(t -> -QuadLL2(Y, X, Z, w1, w2,
                           [t[1],t[2],t[3],t[4:(3+size(X,2))],
                            t[(4+size(X,2)):(3+size(X,2)+size(Z,2))],
                            t[(4+size(X,2)+size(Z,2))]]), θ₀,
             method = BFGS(), f_tol = 1e-10, g_tol = 1e-10).minimizer
