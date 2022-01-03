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
    car_data[!, :Model_id] = Int.(car_data.Model_id)
    # Create unique identifier for each product_marker
    car_data[!, :id] = string.(car_data[!, :Year]) .* "_" .* lpad.( car_data[:, :Model_id], 4, "0")
    # Sort by id
    sort!(car_data, [:id])
    # instruments
    instruments = DataFrame(StatFiles.load(path_to_dir*"/Car_demand_iv_spec1.dta"))
    dropmissing!(instruments)
    instruments[!, :Year] = Int.(instruments.Year)
    instruments[!, :Model_id] = Int.(instruments.Model_id)
    # Create unique identifier for each product_marker
    instruments[!, :id] = string.(instruments[!, :Year]) .* "_" .* lpad.( instruments[:, :Model_id], 4, "0")
    # Sort by id
    sort!(instruments, [:id])
    # simulated income (random coefficient)
    income = DataFrame(StatFiles.load(path_to_dir*"/Simulated_type_distribution.dta"))
    dropmissing!(income)

    return car_data, instruments, income
end

# Function to construct model from data
function construct_model(model_specs::Dict, car_data::DataFrame, instruments::DataFrame, income::DataFrame)
    

    parameters = Primitives()

    # Get id's of the products and markets
    market_id = model_specs[:ids][:market]
    product_id = model_specs[:ids][:product]


    # Create co-variate matrix
    X = car_data[!, model_specs[:covariates]] |> Matrix

    # Create exogenous variables matrix
    Z_exo = car_data[!, model_specs[:exogenous]] |> Matrix

    # Create instruments matrix
    Z_inst = instruments[!, model_specs[:instruments]] |> Matrix

    # Merge exogenous and instruments
    Z = hcat(Z_exo, Z_inst)

    # Create income vector
    Y = income[:, model_specs[:sim_data]] |> Matrix

    # Create a DataFrame for inverse demand estimation 

    inv_dem_est = car_data[!, [:id, :Year, :share, :price]]

    # Initial guess for the inverse demand
    inv_dem = zeros(size(car_data)[1])
    
    inv_dem_est[:, :δ] = inv_dem

    # Initialize GMM Residuals
    gmm_residuals = zeros(size(car_data)[1])
    
    return  Model(parameters, market_id, product_id, X, Z, Y, inv_dem_est, gmm_residuals)
end