#==
    This file conducts the analyses for JF's PS1
==#

using CSV, DataFrames, Optim
#include("./PS1b/JuliaCode/functions.jl")
include("../JuliaCode/functions.jl")
## load the mortgage data as a DataFrame
#df = DataFrame(CSV.File("./PS1b/data/mortgage.csv"))
df = DataFrame(CSV.File("C:/Users/ryana/OneDrive/Documents/School/PhD Economics/Research/GitHub/ECON-899/PS1b/data/mortgage.csv"))

#df[!, :i_25] = df[!, :i_open_year2] .- df[!, :i_open_year5]
#I think i_open_year2-i_open_year5 is Stata notation for
#"include variables i_open_year2 through i_open_year5"
# Yes, it is -Danny

## Separate data into independent variable matrix X and
## dependent variable vector Y
X = df[!,[:i_large_loan,:i_medium_loan,:rate_spread,
          :i_refinance,:age_r,:cltv,:dti, :cu,
          :first_mort_r,:score_0,:score_1, :i_FHA,
          :i_open_year2,:i_open_year3, :i_open_year4,
          :i_open_year5]] |> Matrix;

Y = df[!, :i_close_first_year]; #|> Matrix


## 1. Evaluate functions at β₀ = -1 and β = 0
β = [-1; zeros(size(X, 2), 1)];
LL = likelihood(β, Y, X);
gβ = score(β, Y, X);
H  = Hessian(X, β);

## 2. Compare score and hessian from (1) with numerical
##    first and second derivative of the log-likelihood
gβ_num=∂F(β,Y,X)
    diff_gβ=gβ.-gβ_num
H_num=Find_H_num(β,Y,X)
    diff_H=H-H_num


## 3. Write a routine that solves the maximum likelihood
##    using a Newton algorithm
@time β_Newton = NewtonAlg(Y, X); #Newton(Y, X; β₀ = β);


## 4. Compare the solution and speed with  BFGS and Simplex
f(b) = likelihood(b, Y, X);
@time β_BFGS    = optimize(f, β, BFGS())
@time β_simplex = optimize(f, β, )
