# Load packages
using StatFiles, DataFrames

include("./functions.jl")

# Function to load data
function load_data(path_to_dir::String)
    # Load the data
    # car characteristics
    car_data = DataFrame(StatFiles.load(path_to_dir*"/Car_demand_characteristics_spec1.dta"))
    dropmissing!(car_data)
    car_data[!, :Year] = Int.(car_data.Year)
    # instruments
    instruments = DataFrame(StatFiles.load(path_to_dir*"/Car_demand_iv_spec1.dta"))
    dropmissing!(instruments)
    # simulated income (random coefficient)
    income = DataFrame(StatFiles.load(path_to_dir*"/Simulated_type_distribution.dta"))
    dropmissing!(income)

    return car_data, instruments, income
end

# Function to construct model from data
function construct_model(car_data::DataFrame, instruments::DataFrame, income::DataFrame)
    
    parameters = Primitives()

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

        # Create matrix of product characteristics
        

        # Create zero (or initial guess) inverse demands for the year
        inv_dem = zeros(length(pro))
        demand[year] = Dict(pro .=> inv_dem)

    end

    return  Model(parameters, years, products, shares, prices, income.Var1, X, Z, demand)
end