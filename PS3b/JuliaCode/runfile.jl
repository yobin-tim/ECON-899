#==
    This file conducts the analyses for JF's PS3
==#

# Load the necessary packages
using StatFiles, DataFrames

# Indlude the functions
include("./functions.jl")

# Load the data
# car characteristics
car_data = DataFrame(StatFiles.load("./PS3b/data/Car_demand_characteristics_spec1.dta"))
car_data[!, :Year] = Int.(car_data.Year)
# instruments
instruments = DataFrame(StatFiles.load("./PS3b/data/Car_demand_iv_spec1.dta"))
# simulated income (random coefficient)
income = DataFrame(StatFiles.load("./PS3b/data/Simulated_type_distribution.dta"))
dropmissing!(income)

# Parameters
λₚ = 0.6

# Construct the model
years = unique(car_data.Year)
products = Dict()
shares = Dict()
prices = Dict()
demand = Dict()

for year in years

    # Get the data for the year
    year_data = dropmissing( car_data[car_data.Year .== year, :] )
    pro = year_data.Model_id
    sha = year_data.share
    pri = year_data.price
    products[year] = pro
    shares[year] = Dict(pro .=> sha)
    prices[year] = Dict(pro .=> pri)

    # Create zero (or initial guess) inverse demands for the year
    inv_dem = zeros(length(pro))
    demand[year] = Dict(pro .=> inv_dem)
end


model = Model(λₚ, years, products, shares, prices, income.Var1, demand)

year = 1985

# Get the product indexes
J = model.J[year]
# Get the observed market shares
S = [model.S[year][j] for j in J]
# Get the observed prices
P = [model.P[year][j] for j in J]
# Get the income levels
Y = model.Y
# Get the inital guess for the inverse demand
δ = [model.δ[year][j] for j in J]

μ = repeat(Y', length(J), 1) .* P


choice_probability(δ, μ)

inverse_demand_cm(δ, S, μ, 0.01; max_iter = 1000)