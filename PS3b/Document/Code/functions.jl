# Load packages
using LinearAlgebra, Parameters, Optim

# Auxiliary functions
# Choice probability function
function choice_probability(δ::Array{Float64}, μ::Array{Float64}; eval_jacobian::Bool = false)

    # number of individuals and choicesm
    J, R = size(μ)

    # Compute choice probabilities
    Λ = exp.(δ .+ μ)
    Σ = Λ ./ (1 .+ sum(Λ, dims = 1))
    σ = sum(Σ, dims = 2) / R

    if eval_jacobian
        # Compute Jacobian
        Δ = (1 / R) * ((I(J) .* (Σ * (1 .- Σ)')) - ((1 .- I(J)) .* (Σ * Σ'))) ./ σ
        return σ, Δ
    else
        return σ, nothing
    end

end

# segment_data by market for demand estimation
function segment_data(model, market)

    # Get market id column
    market_id_col = model.market_id

    # Filter data by market
    data = model.inv_dem_est[model.inv_dem_est[!, market_id_col] .== market, :]

    # Get the observed market shares
    S = data.share
    # Get the observed prices
    P = data.price
    # Get the income levels
    Y = model.Y
    # Get the inital guess for the inverse demand
    δ = data.δ

    return S, P, Y, δ
end

# Model Structures
# Primitives
@with_kw struct Primitives
    λₚ_range :: Array{Float64} = [0, 1]
end


# Model
mutable struct Model
    # Parameters
    parameters      ::Primitives                # Parameters of the model

    # Data
    market_id       :: Any                      # Market id column
    product_id      :: Any                      # Product id column
    X               :: Array{Float64, 2}        # Matrix of covariates
    Z               :: Array{Float64, 2}        # Matrix of instruments
    Y               :: Array{Float64, 2}        # Matrix of simulated data
    inv_dem_est     :: DataFrame                # DataFrame of for demand estimation

    # GMM estimation
    ρ               :: Array{Float64}           # GMM Residuals
end



# Demand inverter
function inverse_demand(model::Model, λₚ::Float64, market; method::String="Newton", max_iter = Inf, verbose::Bool = false, recalc_δ::Bool = true)
    if method=="Newton" #Otherwise this prints too often when we do the contraction mapping
        print("\n")
        print("Market: $(market)")
        print("\n")
    end
    # Check the method
    valid_methods = ["Newton", "Contraction Mapping"]
    @assert (method ∈ valid_methods)

    # Get the data
    S, P, Y, δ = segment_data(model, market)

    if recalc_δ
        # Recalculate the initial guess for the inverse demand
        δ = zeros(size(δ))
    end

    # Compute the matrix μ[i,j] = λₚ * Y[i] * P[j]
    μ = λₚ * repeat(Y', length(S), 1) .* P

    # Compute the inverse demand

    # Initial guess for the inverse demand
    δ₀ = copy( δ )
    δ₁ = copy( δ )

    # Initialize the iteration counter
    iter = 0

    # Initialize the error
    err = 100
    eval_jacobian = false

    ε = 1e-12
    ε₁ = ( method == "Newton" ) ? 1 : -Inf

    # Iterate until convergence

    err_list = []
    method_flag = "Contraction Mapping"

    while (err > ε) && (iter < max_iter)
        # Compute the choice probability and the Jacobian
        if (method == "Newton") && (err < ε₁)
            eval_jacobian = true
            method_flag = "Newton"
        else # This will bring it back to contraction mapping if it diverges
            eval_jacobian = false
            method_flag = "Contraction Mapping"
        end

        σ, Δ = choice_probability(δ₀, μ, eval_jacobian=eval_jacobian)

        # Compute the inverse demand
        if (method == "Newton") && (err < ε₁)
            #I added the ./S after talking with Michael Nattinger
            #It also lines up with JF's ox code in blp_func_ps.ox
            #(From Danny) JF divides by "Shat", or predicted shares,
            # aka σ
            #δ₁ = δ₀ + inv(Δ ./ S) * (log.(S) - log.(σ))
            δ₁ = δ₀ + inv(Δ) * (log.(S) - log.(σ))
        else
            δ₁ = δ₀ + log.(S) - log.(σ)
        end

        # Update the error
        #err = maximum( abs.(δ₁ - δ₀) )
        err = norm(δ₁ - δ₀)

        push!(err_list, err)
        # Update the inverse demand
        δ₀ = copy(δ₁)

        # Update the iteration counter
        iter = iter + 1
        if (iter % 10 == 0) && verbose
            println("Iteration = $iter, Method = $method_flag , error = $err, tolerance = $ε, error > tolerance = $(err > ε)")
        end

    end
    # println("Iteration = $iter, Method = $method_flag, error = $err, tolerance = $ε, error > tolerance = $(err > ε), θ = $θ")
    market_id_col = model.market_id
    model.inv_dem_est[model.inv_dem_est[!, market_id_col] .== market, :δ] .= δ₁[:, 1]

    return err_list

end


function gmm(model, λ; Return_ρ::Bool=false, SpecifyW::Bool=false,
        SpecifiedW::Array{Float64, 2}=zeros(2,2))
    markets = unique(model.inv_dem_est[!, model.market_id])
    for  market in markets
        inverse_demand(model, λ, market, method = "Contraction Mapping")
    end

    # Iv regression
    X = model.X
    Z = model.Z
    W = inv(Z'Z)
    if SpecifyW
        W=SpecifiedW
    end
    δ = model.inv_dem_est.δ
    β_iv = inv((X'Z)*W*(Z'X))*(X'Z)*W*(Z'δ)

    ρ = (δ - X*β_iv)
    if ~Return_ρ
        return ρ'Z*W*Z'*ρ/100
    else
        return ρ/100
    end
end


function TwoStage_gmm(model)
    λhat = optimize(λ -> gmm(model, λ[1]), [.6],
                 method = BFGS(), f_tol = 1e-5, g_tol = 1e-5).minimizer
    ξhat=gmm(model, λhat[1],Return_ρ=true)
    #OptimalW=inv( (model.Z * ξhat)*transpose(model.Z * ξhat) )
    #Maybe it ought to be this?
    OptimalW=inv( (model.Z .* ξhat)'*(model.Z .* ξhat) )
    #OptimalW=inv( dot(model.Z, ξhat)transpose(dot(model.Z, ξhat) ))
    #print(OptimalW)
    print(λhat)
    λhat_SecondStage=optimize(λ -> gmm(model, λ[1],SpecifyW=true,SpecifiedW=OptimalW),
                λhat, method = BFGS(), f_tol = 1e-5, g_tol = 1e-5).minimizer
    return λhat_SecondStage
end
