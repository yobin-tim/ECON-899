# Load packages
using LinearAlgebra, Parameters, Optim, CSV, Printf

## Structure that contains the parameters of the model
@with_kw mutable struct Primitives
    λ₀::Int64   = -4
    α::Int64    = 2
    β::Float64  = 0.99
    pᵣ::Int64   = 4
    pₛ::Int64   = 1
    imax::Int64 = 8
    λ::Float64 

    # Markov transition matrix
    Π::Array{Float64,2} = [
        0.9 0.1
        0.9 0.1
    ]

    # Utility function
    U::Function = (a, i, c, p, λ, ϵ) -> (a == 1) ? α * c - p + ϵ : (
                        (i > 0) ? α * c + ϵ[1] : λ * (c > 0) + ϵ)

    # Parameters that must be loaded:
    # State space S 
    # Transition matrix F(s'|s, a = 0)
    # Transition matrix F(s'|s, a = 1)
    S::Array{Int64,2} = (DataFrame(CSV.File(
        "./PS4b/ox_code/PS4_state_space.csv"))|>Matrix)[:, 3:end]
    F₀::Array{Float64,2} = (DataFrame(CSV.File(
        "./PS4b/ox_code/PS4_transition_a0.csv"))|>Matrix)[:, 3:end]
    F₁::Array{Float64,2} = (DataFrame(CSV.File(
        "./PS4b/ox_code/PS4_transition_a1.csv"))|>Matrix)[:, 3:end]
end # Primitives struct 


## Function calculates expected value
function ExpVal(prim::Primitives, P, θ; primβ = true)

    # unpack relevant primitives
    @unpack α, β, S, U, F₀, F₁ = prim
    γ = MathConstants.eulergamma;

    # replace β with second argument of θ is primβ not specified
    if (primβ); λ = θ; else; λ, β = θ[1], θ[2]; end;

    # define utility levels for each state space and choice of a 
    U₁ = U.(1, S[:, 1], S[:, 2], S[:, 3], λ, 0)
    U₀ = U.(0, S[:, 1], S[:, 2], S[:, 3], λ, 0)

    # steps from JF's code (not totally clear to me)
    ϵ₀, ϵ₁ = γ .- log.(1 .- P), γ .- log.(P)
    F = F₀.*(1 .- P) .+ F₁.*P

    EU = (1 .- P).*(U₀ .+ ϵ₀) .+ P.*(U₁ .+ ϵ₁)
    EV = inv(I(size(P, 1)) .- β.*F)*EU

    V = [(U₀ + β*F₀*EV) (U₁ + β*F₁*EV)]

    return EV, exp.(V[:, 2]) ./ sum(exp.(V), dims = 2)

end # ExpVal()

## Wrapper function for the expected value function that performs
## the CCP algorithm
function CCP(prim::Primitives, P₀; ε = 10e-10, err = 100, N = 100, trueλ = true)

    # unpack relevant primitives
    if (trueλ)
        @unpack λ₀ = prim
    else
        @unpack λ₀, λ = prim; λ₀ = λ
    end


    # initialize iteration counter
    i = 1;

    # iterate the expected value function until convergence
    while err > ε
        global EV, P = ExpVal(prim, P₀, λ₀)
        err = norm(P - P₀)

        # print progress every N iterations
        if (i % N == 0)
            println(@sprintf "%.0f iterations, err = %.5f" i err)
        end
        i += 1; P₀ = P;
    end # convergence loop

    # return converged choice probability and expected value
    return EV, P

end # CCP()

## Function that evaluates the log-likelihood for choices aᵢ at a given λ and β
function LL(aᵢ, Sᵢ, θ, P₀; primβ = true)

    # use the CCP to solve for P given θ
    if (primβ); prim.λ = θ; else; prim.λ, prim.β = θ[1], θ[2]; end;

    # solve for fixed-point CCP
    P = CCP(prim, P₀; trueλ = false)[2]

    # calculate likelihood for each observation
    Lᵢ = zeros(length(unique(Sᵢ)), 1)
    for i ∈ 0:(length(unique(Sᵢ)) - 1)
        Lᵢ[i+1] = sum(aᵢ[Sᵢ .== i].*log.(P[i + 1]) + 
                    (1 .- aᵢ[Sᵢ .== i]).*log.(1 .- P[i + 1]))
    end

    # return log-likelihood
    return sum(Lᵢ)
end # log-likelihood function, LL()
