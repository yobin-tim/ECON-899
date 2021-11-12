#==
    This file conducts the analyses for JF's PS2
==#

using StatFiles, DataFrames, Optim, BenchmarkTools, Latexify, CSV
# We can use BenchmarkTools for better precision. Just need to replace time
# with btime. The runtime of the overall code get longer as btime runs the
# code multiple times to reduce noise

include("./PS2b/JuliaCode/functions.jl")
#include("./functions.jl")

## load the mortgage data and sparse grid weights as a DataFrames (and 
## convert weights to matrices
df = DataFrame(StatFiles.load("PS2b/data/Mortgage_performance_data.dta"))  
w1 = DataFrame(CSV.File("PS2b/data/KPU_d1_l20.csv")) |> Matrix
w2 = DataFrame(CSV.File("PS2b/data/KPU_d2_l20.csv")) |> Matrix

# Use this if you are loading data from the root folder.
#df = DataFrame(CSV.File("../data/mortgage_performance_data.csv"))
#df = DataFrame(CSV.File("C:/Users/ryana/OneDrive/Documents/School/PhD Economics/Research/GitHub/ECON-899/PS1b/data/mortgage.csv"))

## Separate data into independent variable matrices X and Z and
## dependent variable vector Y
X = df[!, [:score_0, :rate_spread, :i_large_loan, :i_medium_loan,
            :i_refinance, :age_r, :cltv, :dti, :cu,
            :first_mort_r, :i_FHA, :i_open_year2,
            :i_open_year3, :i_open_year4, :i_open_year5]] |> Matrix;

Z = df[!, [:score_0, :score_1, :score_2]] |> Matrix;

Y = df[!, :duration]; #|> Matrix


<<<<<<< HEAD
for name in names(df)
    println(name)
end
## 1. Evaluate log-likelihood using the quadrature method
println("See likelihood() function")

## 2. Evaluate simulated log-likelihood function using GHK
θ = [0, 0, 0, zeros(size(X, 2), 1), zeros(size(Z, 2), 1), .5];
@btime β_BFGS = optimize(t->-likelihood(Y, X, Z, w1, w2, t),
                            θ, method=BFGS(), 
                            f_tol=1e-32,g_tol=1e-32).minimizer


## 3. Evaluate simulated log-likelihood function using accept/reject


## 4. Compare predicted choice probabilities for each of the above methods
θ₀ = [0, -1, -1, 0*ones(size(X, 2), 1), 0.3*ones(size(Z, 2), 1), 0.5]

likelihood(Y, X, Z, w1, w2, θ₀)