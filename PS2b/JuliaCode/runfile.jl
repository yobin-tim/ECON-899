#==
    This file conducts the analyses for JF's PS2
==#

using CSV, DataFrames, Optim, BenchmarkTools, Latexify
# We can use BenchmarkTools for better precision. Just need to replace time
# with btime. The runtime of the overall code get longer as btime runs the
# code multiple times to reduce noise

# include("./PS1b/JuliaCode/functions.jl")
include("./functions.jl")

## load the mortgage data as a DataFrame
#df = DataFrame(CSV.File("../data/mortgage_performance_data.csv"))

# Use this if you are loading data from the root folder.
df = DataFrame(CSV.File("PS2b/Data/mortgage_performance_data.csv"))  
#df = DataFrame(CSV.File("C:/Users/ryana/OneDrive/Documents/School/PhD Economics/Research/GitHub/ECON-899/PS1b/data/mortgage.csv"))

## Separate data into independent variable matrixes X and Z and
## dependent variable vector Y
X = df[!, [:score_0, :rate_spread, :i_large_loan, :i_medium_loan,
            :i_refinance, :age_r, :cltv, :dti, :cu,
            :first_mort_r, :i_FHA, :i_open_year2,
            :i_open_year3, :i_open_year4, :i_open_year5]] |> Matrix;

Z = df[!, [:score_0, :score_1, :score_2]] |> Matrix;

Y = df[!, :duration]; #|> Matrix


## 2. Evaluate simulated log-likelihood function using GHK
@btime β_BFGS = optimize(b->-likelihood(b, Y, X), β, method=BFGS(),
                            f_tol=1e-32,g_tol=1e-32).minimizer


## 3. Evaluate simulated log-likelihood function using accept/reject


## 4. Compare predicted choice probabilities for each of the above methods


## 5.