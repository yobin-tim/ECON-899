#==
    This file conducts the analyses for JF's PS1
==#

using CSV, DataFrames, Optim, BenchmarkTools, Latexify
# We can use BenchmarkTools for better precision. Just need to replace time
# with btime. The runtime of the overall code get longer as btime runs the
# code multiple times to reduce noise

# include("./PS1b/JuliaCode/functions.jl")
include("../JuliaCode/functions.jl")
## load the mortgage data as a DataFrame
# df = DataFrame(CSV.File("../data/mortgage.csv"))

# Use this if you are loading data from the root folder.
df = DataFrame(CSV.File("PS1b/Data/mortgage.csv"))  
#df = DataFrame(CSV.File("C:/Users/ryana/OneDrive/Documents/School/PhD Economics/Research/GitHub/ECON-899/PS1b/data/mortgage.csv"))

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


gβ = score(β, Y, X)
# The transpose of the score evaluated at β is 
latexify(gβ')

H  = Hessian(X, β)
# The Hessian evaluated at β is
latexify(H)

## 2. Compare score and hessian from (1) with numerical
##    first and second derivative of the log-likelihood
gβ_num=∂F(β,Y,X)
# gβ_num=score_num(β,Y,X)
diff_gβ=gβ.-gβ_num

H_num=Find_H_num(β,Y,X)
diff_H=H-H_num


## 3. Write a routine that solves the maximum likelihood
##    using a Newton algorithm
@btime β_Newton = NewtonAlg(Y, X); #Newton(Y, X; β₀ = β);
# β_Newton = NewtonAlg(Y, X);

## 4. Compare the solution and speed with  BFGS and Simplex
#f(b) = likelihood(b, Y, X);
#Optimize minimizes the function, so we need to use the negative of liklihood to maximize
println("\n For Quasi-Newton Methods:")
print("\n The BFGS algorithm takes")
# @btime β_BFGS    = optimize(b->-likelihood(b, Y, X), β, BFGS(),abs_tol=1e-12).minimizer
@btime    β_BFGS    = optimize(b->-likelihood(b, Y, X), β, method=BFGS(),
        f_tol=1e-32,g_tol=1e-32).minimizer

print("\n The Simplex algorithm takes")
@btime β_simplex = optimize(b->-likelihood(b, Y, X), β, NelderMead()).minimizer;
    # β_simplex = optimize(b->-likelihood(b, Y, X), β, method=NelderMead(),
    # f_tol=1e-32,g_tol=1e-32).minimizer
