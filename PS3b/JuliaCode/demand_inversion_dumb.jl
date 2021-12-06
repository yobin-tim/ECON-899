using LinearAlgebra
model = construct_model(model_specs, car_data, instruments, income)
S, P, Y, δ = segment_data(model, market)

function ch_prob(δ, P, Y)
    σ = zeros( length(δ), length(Y) )
    σj = zeros( size(δ) )
    for j in 1:length(δ)
        for r in 1:length(Y)
            μ = Y[r] * P[j]
            num = exp(δ[j] + μ)
            den = 1
            for j_prime in 1:length(δ)
                μ = Y[r] * P[j_prime]
                den += exp(δ[j_prime] + μ)
            end
            σj[j] += num / den
            σ[j, r] = num / den
        end 
        σ[j] /= length(Y)
    end
    return σ, σj
end

δ₀ = copy( δ )
δ₁ = copy( δ )

err, iter,err_list = 100, 0, []
while (err > 1) && (iter < 500)
    _, σj = ch_prob(δ₀, P, Y)

    δ₁ = δ₀ + log.(S) - log.(σj)

    # Update the error
    err = maximum( abs.(δ₁ - δ₀) )
    push!(err_list, err)
    # Update the inverse demand
    δ₀ = copy(δ₁)
    
    # Update the iteration counter
    iter = iter + 1

end


# Compute Jacobian
σ, σj = ch_prob(δ₀, P, Y)

R = size(Y)[1]
r = size(δ)[1]

Δ =  1/R *( I(r) .* (σ * (1 .- σ)') - ((1 .- I(r)) .* (σ * σ'))) ./ σj'

inv(Δ) 

iter = 0
while (err > 1e-12) && (iter < 20)
    σ, σj = ch_prob(δ₀, P, Y)

    Δ =  1/R *( I(r) .* (σ * (1 .- σ)') - ((1 .- I(r)) .* (σ * σ'))) ./ σj'
    
    δ₁ = δ₀ - inv(Δ) * (log.(S) - log.(σj))

    # Update the error
    err = maximum( abs.(δ₁ - δ₀) )
    push!(err_list, err)
    # Update the inverse demand
    δ₀ = copy(δ₁)
    
    # Update the iteration counter
    iter = iter + 1

    println(err)
end

σ = ch_prob(δ₀, P, Y)
Δ = 1/R * ((I(r) .* (σ * (1 .- σ)')) - ((1 .- I(r)) .* (σ * σ'))) ./ σ 

δ₁ = δ₀ - inv(Δ) * (log.(S) - log.(σ))

err = maximum( abs.(δ₁ - δ₀) )

δ₀ = copy(δ₁)