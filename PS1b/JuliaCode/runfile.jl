#==
    This file conducts the analyses for JF's PS1
==#

using CSV, DataFrames

## load the mortgage data as a DataFrame 
df = DataFrame(CSV.File("../data/mortgage.csv"))

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
          :i_open_year5]] |> Matrix

Y = df[!, :i_close_first_year] #|> Matrix

