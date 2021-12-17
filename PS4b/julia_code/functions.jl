# Load packages
using LinearAlgebra, Parameters, Optim, CSV

## Structure that contains the parameters of the model
@with_kw mutable struct Primitives
    λ::Int64 = -4
    α::Int64 = 2
    β::Float64 = 0.99
    pᵣ::Int64 = 4
    pₛ::Int64 = 1
    imax::Int64 = 8

    # Markov transition matrix
    Π::Array{Float64, 2} = [
        0.9 0.1
        0.9 0.1
    ]

    # Parameters that must be loaded:
    # State space S 
    # Transition matrix F(s'|s, a = 0)
    # Transition matrix F(s'|s, a = 1)
    S::Array{Int64, 2} = (DataFrame(CSV.File(
        "./PS4b/ox_code/PS4_state_space.csv"))|>Matrix)[:, 3:end]
    F₀::Array{Float64, 2} = (DataFrame(CSV.File(
        "./PS4b/ox_code/PS4_transition_a0.csv")) |> Matrix)[:, 3:end]
    F₁::Array{Float64, 2} = (DataFrame(CSV.File(
        "./PS4b/ox_code/PS4_transition_a1.csv")) |> Matrix)[:, 3:end]
end # Primitives struct 
