# Load packages
using Parameters

# Auxiliary functions
# Choice probability function
function choice_probability(δ::Array{Float64}, μ::Array{Float64}; eval_jacobian::Bool = False)
    
    R = length(μ)
    Δ = exp.( repeat(δ, (1, lenght(μ))) .* repeat(μ', (lenght(δ), 1)) )
    Δ = Δ./(1 + sum(Δ, dims=2))

    if eval_jacobian
        #return Δ, jacobian(Δ, δ)
    else
        return sum(Δ, dims=1)/R
    end
end

# Model Structure
mutable struct Model 

    # Parameters
    λₚ :: Float64 # Somehting

    # Data
    
    # Functions
    # choice_probability :: Function(Float64, Float64) -> Float64

    # Demand

end




# Demand inversion Contraction Mapping (Barry et al., 1995)
function inverse_demand_cm(δ, s, μ, tol; max_iter::Int64 = 1000)

    # Initial guess for the inverse demand
    δ₀ = copy(δ)

    # Initialize the iteration counter
    iter = 0
    # Initialize the error
    err = 100

    # Iterate until convergence
    while (err < tol) && (iter < max_iter)
        # Compute the choice probability
        σ = choice_probability(δ₀, μ)

        # Compute the inverse demand
        δ = δ₀ + log.(s) - log.(σ)

        # Update the inverse demand
        δ₀ = copy(δ)

        # Update the error
        err = norm(δ₁ - δ₀)

        # Update the iteration counter
        iter = iter + 1
    end

    return δ
end

# Demand inverter
function inverse_demand(data, year, params::Array{Float64} ; method::String="Newton")

    if method=="Newton"

    elseif method=="CM"

    else
        error("Method not implemented, implemented methods are \"Newton\" and \"CM\" ")
    end

end


inverse_demand([1.0]; method="CM")