#==
This file defines functions used in JF's PS2
==#
using Optim, Distributions, Parameters, LinearAlgebra

# structure of model parameters
mutable struct ModelParameters
    α₀::Float64
    α₁::Float64
    α₂::Float64
    β::Array{Float64}
    γ::Array{Float64}
    ρ::Float64
end # parameters struct

# Calculate log-likelihood using quadrature method
function QuadLL2(Y, X, Z, W1, W2, θ)

    u = W1[:, 1]; w = W1[:, 2]
    μ₀ = W2[:, 1]; μ₁ = W2[:, 2]; ω = W2[:, 3]

    param = ModelParameters(θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
    @unpack α₀, α₁, α₂, β, γ, ρ = param

    # Calculate σ₀ and σ₀²
    σ₀² = 1/(1-ρ)^2
    σ₀  = 1/(1-ρ)

    tmp = α₀ .+ X*β .+ Z*γ;
    mρ = zeros(size(X,1), size(u,1));

    for i in 1:size(X,1) # For each observation, get range based on domain (0,1) at t = 0
        mρ[i,:] = log.(u') .+ tmp[i]
    end

    tmp = α₁ .+ X*β .+ Z*γ
    mρ1 = zeros(size(X,1), size(u,1));

    for i in 1:size(X,1) # For each observation, get range based on domain (0,1) at t = 1
        mρ1[i,:] = log.(u') .+ tmp[i]
    end

    mdρ = ones(size(mρ,1), size(mρ,2));

    for i in 1:size(X,1) # For each observation, get Jacobian 1/u
        mdρ[i,:] = mdρ[i,:] ./u
    end

    L1 = cdf.(Normal(), (-α₀ .- X*β .- Z*γ)./σ₀)

    density = pdf.(Normal(), mρ./σ₀)./σ₀

    L2 = (cdf.(Normal(), - α₁ .- X*β .- Z*γ .- ρ .* mρ) .* density.* mdρ) * w

    density = pdf.(Normal(), mρ1 - ρ*mρ) .* pdf.(Normal(), mρ./σ₀) ./ σ₀

    L3 = (cdf.(Normal(), - α₂ .- X*β .- Z*γ .- ρ .* mρ1) .* density .* mdρ .* mdρ) * w

    L4 = 1 .- L1 .- L2 .- L3

    ll = 0

    for i = 1:size(Y, 1)

        if Y[i] == 1

            ## If the likelihood becomes minus, I evaluate this value as 1e-10.
            if L1[i] < 0
                L1[i]=1e-10
            else
            end

            ll = ll + log(L1[i])

        elseif Y[i] == 2

            if L2[i] < 0
                L2[i]=1e-10
            else

            end

            ll = ll + log(L2[i])

        elseif Y[i] == 3

            if L3[i] < 0
                L3[i]=1e-10
            else
            end

            ll = ll + log(L3[i])

        elseif Y[i] == 4

            if L4[i] < 0
                L4[i]=1e-10
            else
            end

            ll = ll + log(L4[i])
        end

    end # for i

    return(ll)

end

# Calculate log-likelihood using quadrature method
function QuadLL(Y, X, Z, W1, W2, θ)
    # separate weights and nodes from W1 and W2
    u = W1[:, 1]; w = W1[:, 2]
    μ₀ = W2[:, 1]; μ₁ = W2[:, 2]; ω = W2[:, 3]

    # unpack model parameters
    param = ModelParameters(θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
    @unpack α₀, α₁, α₂, β, γ, ρ = param

    # Calculate σ₀ and σ₀²
    σ₀² = 1/(1-ρ)^2
    σ₀  = 1/(1-ρ)

    # map integral bounds to [0, 1] Upper bound (-Infinity, α + Xβ + Zγ)
    # ρ(u) = ln(u) + α + Xβ + Zγ

    b₀ = (a, x, z) -> log.(a) .+ (α₀ + dot(x, β) + dot(z, γ))
    b₁ = (a, x, z) -> log.(a) .+ (α₁ + dot(x, β) + dot(z, γ))

    # per-observation likelihood:

    L1 = (x, z) -> log(cdf(Normal(), (-α₀ - dot(x, β) - dot(z, γ))/σ₀))
    # I think this should just be σ₀ not σ₀²

    # I think the location of the () about σ₀ wrong.
    #L2 = (x, z) -> log(sum(w.*
    #    (cdf.(Normal(), (-ρ)*b₀(u, x, z) .- (α₁ .+ dot(x, β) + dot(z, γ)))./σ₀).*
    #    pdf.(Normal(), b₀(u./σ₀, x, z)) ./ u))

    function L2(x, z)
        out=0
        try
            out=log(sum(w.*
                (cdf.(Normal(), -α₁ .- dot(x, β) .- dot(z, γ) .-ρ*b₀(u, x, z))).* # Function
                (pdf.(Normal(), b₀(σ₀ .* u, x, z)./σ₀)./σ₀) .* # Density
                (1 ./u))) # Jacobian

        catch #Sometimes parameters will be tried that make the above try to take the log
            #of a negative number. This is an attempted fix for that which just returns something
            #awful in that case so that it won't be picked
            out=log(1e-5)
            print("Attempted a point that does not work.")
        end
        return out
    end


    function L3(x, z)
        out=0
        try
            out=log(sum(ω.*((cdf.(Normal(), (-ρ)*b₁(μ₁, x, z) .- (α₂ .+ dot(x, β) +
                dot(z, γ))))./σ₀).*pdf.(Normal(), b₀(μ₀, x, z)./σ₀).*pdf.(Normal(), b₁(μ₁, x, z).-
                ρ*b₀(μ₀, x, z)) ./ (μ₀ .* μ₁)))
        catch
            out=log(1e-5)
            print("Attempted a point that does not work.")
        end
        return out
    end
    function L4(x, z)
        out=0
        try
            out=log(sum(ω.*(cdf.(Normal(), -(ρ)*b₁(μ₁, x, z) .+ (α₂ + dot(x, β) +
                dot(z, γ)))./σ₀).*pdf.(Normal(), b₁(μ₁./σ₀, x, z)).*pdf.(Normal(), b₁(μ₁, x, z) .-
                ρ*b₀(μ₀, x, z)) ./ (μ₀ .* μ₁)))
        catch
            out=log(1e-5)
            print("Attempted a point that does not work.")
        end
        return out
    end

    # calculate the log-likelihood for all observations
    ll = 0
    for i = 1:size(Y, 1)
        if Y[i] == 1
            ll = ll + L1(X[i, :], Z[i, :])
        elseif Y[i] == 2
            ll = ll + L2(X[i, :], Z[i, :])
        elseif Y[i] == 3
            ll = ll + L3(X[i, :], Z[i, :])
        elseif Y[i] == 4
            ll = ll + L4(X[i, :], Z[i, :])
        end
    end # for i

    return ll
end # quadrature log-likelihood function



# Calculate log-likelihood using quadrature method
function GHKLL(Y, X, Z, θ; sims = 100)

    # unpack model parameters
    param = ModelParameters(θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
    @unpack α₀, α₁, α₂, β, γ, ρ = param
    σ₀  = 1/(1-ρ)

    ll = 0
    for i=1:size(Y, 1)
        ll_i = 1
        ϵ_draws = zeros(sims, 3)
        if Y[i] > 1 # Need to draw from a distribution which won't make the borrower repay in period 1
            ϵ_draws[:, 1] = rand.(truncated(Normal(0, σ₀), -Inf, -α₀ - dot(X[i, :], β) - dot(Z[i, :], γ)), sims)
        elseif Y[i] == 1 #Find the probability that this draw would have occured
            ll_i=1-cdf(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ))/σ₀)
        end
        if Y[i] > 2 # Need to draw from a distribution which won't make the borrower repay in period 2
            ϵ_draws[:, 2] = [rand(truncated(Normal(0, σ₀), -Inf, -α₀ - dot(X[i, :], β) - dot(Z[i, :], γ) - ρ*ϵ_draws[si, 1])) for si = 1:sims]
        elseif Y[i] == 2 # Find the probability that this draw would have occured
            ll_i = (1/sims)*cdf(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ))/σ₀)*
                sum(1 .- cdf.(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ)) .- ρ*ϵ_draws[:, 1]))
        end
        if Y[i] > 3 # Need to draw from a distribution which won't make the borrower repay in period 3
            ϵ_draws[:, 3]=[rand(truncated(Normal(0, σ₀), -Inf, -α₀ - dot(X[i, :], β) - dot(Z[i, :], γ) - (ρ^(2))*ϵ_draws[si, 1] - ρ*ϵ_draws[si, 2])) for si = 1:sims]
            # Find the probability that Y[i]==4 would have occured
            ll_i = (1/sims)*cdf(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ))/σ₀)*
                sum((cdf.(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ)) .- ρ*ϵ_draws[:, 1])).*
                cdf.(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ)) .- (ρ^(2))*ϵ_draws[:, 1] .- ρ*ϵ_draws[:, 2]))
        elseif Y[i] == 3 # Find the probability that this draw would have occured
            ll_i = (1/sims)*cdf(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ))/σ₀)*
                sum(  (cdf.(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ)) .- ρ*ϵ_draws[:, 1])).*
                (1 .- cdf.(Normal(), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ)) .- (ρ^(2))*ϵ_draws[:, 1] .- ρ*ϵ_draws[:, 2]))   )
        end
        ll += log(ll_i)
    end
    return ll
end # GHK log-likelihood function

function GHKLL2(Y, X, Z, θ; sims = 100)
    #ref: http://fmwww.bc.edu/repec/bocode/g/GHK_note.pdf
    param = ModelParameters(θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
    @unpack α₀, α₁, α₂, β, γ, ρ = param
    σ₀  = 1/(1-ρ)

    ll_store = 0

    for i=1:size(Y, 1)

        if Y[i] == 1
            ll = cdf(Normal(0, 1), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ))/σ₀)
        elseif Y[i] == 2
            η₀ = rand.(truncated(Normal(0, 1), -Inf, (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./σ₀), 100)
            ll = mean(cdf(Normal(0, 1), (-α₁ - dot(X[i, :], β) - dot(Z[i, :], γ) .- η₀)./(ρ .* σ₀)) .*
                cdf(Normal(0, 1), (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./(σ₀)))
        elseif Y[i] == 3
            η₀ = rand.(truncated(Normal(0, 1), -Inf, (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./σ₀), 100)
            tmp = ((α₁ + dot(X[i, :], β) + dot(Z[i, :], γ) .- η₀)./σ₀)
            η₁ = zeros(100)
            for i in 1:100
                η₁[i] = (rand.(truncated(Normal(0, 1), -Inf, (tmp[i])[1]), 1))[1]
            end
            ll = mean(cdf(Normal(0, 1), (-α₂ - dot(X[1, :], β) - dot(Z[i, :], γ) .- η₀ .- ρ.* η₁)./((ρ.^2) .* σ₀)) .*
                cdf(Normal(0, 1), (α₁ + dot(X[i, :], β) + dot(Z[i, :], γ) .- η₁)./(ρ .* σ₀)) .*
                cdf(Normal(0, 1), (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./(σ₀)))
        else

            ll1 = cdf(Normal(0, 1), (-α₀ - dot(X[i, :], β) - dot(Z[i, :], γ))/σ₀)

            η₀ = rand.(truncated(Normal(0, 1), -Inf, (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./σ₀), 100)
            ll2 = mean(cdf(Normal(0, 1), (-α₁ - dot(X[i, :], β) - dot(Z[i, :], γ) .- η₀)./(ρ .* σ₀)) .*
                cdf(Normal(0, 1), (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./(σ₀)))

            η₀ = rand.(truncated(Normal(0, 1), -Inf, (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./σ₀), 100)
            tmp = ((α₁ + dot(X[i, :], β) + dot(Z[i, :], γ) .- η₀)./σ₀)
            η₁ = zeros(100)
            for j in 1:100
                η₁[j] = (rand.(truncated(Normal(0, 1), -Inf, (tmp[j])[1]), 1))[1]
            end
            ll3 = mean(cdf(Normal(0, 1), (-α₂ - dot(X[i, :], β) - dot(Z[i, :], γ) .- η₀ .- ρ.* η₁)./((ρ.^2) .* σ₀)) .*
                cdf(Normal(0, 1), (α₁ + dot(X[i, :], β) + dot(Z[i, :], γ) .- η₁)./(ρ .* σ₀)) .*
                cdf(Normal(0, 1), (α₀ + dot(X[i, :], β) + dot(Z[i, :], γ))./(σ₀)))

            ll = 1 - ll1 - ll2 - ll3

        end

        ll_store += log.(ll)

    end

    return(ll_store)

end



# Calculate log-likelihood using accept-reject method
function AcceptRejectLL(Y, X, Z, θ; sims = 100, k = maximum(Y))

    # unpack model parameters
    param = ModelParameters(θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
    @unpack α₀, α₁, α₂, β, γ, ρ = param
    σ₀ = 1 / (1 - ρ)

    ll = 0 # initialize log-likelihood

    # Define index functions for each outcome
    I1 = (x, z, ε) -> ε .< -(α₀ .+ x * β .+ z * γ)
    I2 = (x, z, ε) -> (ε[:, 1] .< -(α₀ .+ x * β .+ z * γ)) .& (ε[:, 2] .< -(α₁ .+ x * β .+ z * γ) .- ρ * ε[:, 1])
    I3 = (x, z, ε) -> (ε[:, 1] .< -(α₀ .+ x * β .+ z * γ)) .& (ε[:, 2] .< -(α₁ .+ x * β .+ z * γ) .- ρ * ε[:, 1]) .& (
        ε[:, 3] .< -(α₁ .+ x * β .+ z * γ) .- (ρ^2) * ε[:, 1] .- ρ * ε[:, 2])
    I4 = (x, z, ε) -> (ε[:, 1] .< -(α₀ .+ x * β .+ z * γ)) .& (ε[:, 2] .< -(α₁ .+ x * β .+ z * γ) .- ρ * ε[:, 1]) .& (
        ε[:, 3] .< α₁ .+ x * β .+ z * γ .- (ρ^2) * ε[:, 1] .- ρ * ε[:, 2])

    # Calculate log-likelihood for Y = 1 observations
    x, z = repeat(X[Y.==1, :], inner = [sims, 1]), repeat(Z[Y.==1, :], inner = [sims, 1])
    ε = rand.(Normal(0, σ₀), size(x, 1))
    for i = 1:sum(Y .== 1)
        ind = ((i-1)*sims+1):(i*sims)
        if sum(I1(x[ind, :], z[ind, :], ε[ind, :]))!=0
            ll += log(1/sum(I1(x[ind, :], z[ind, :], ε[ind, :])))
            #ll += log(sum(I1(x[ind, :], z[ind, :], ε[ind, :])) / sims)
        else
            ll+=log(1/sims) #Act as though it happened at least once to avoid -infinity
        end
    end

    # Calculate log-likelihood for Y = 2 observations
    x, z = repeat(X[Y.==2, :], inner = [sims, 1]), repeat(Z[Y.==2, :], inner = [sims, 1])
    ε = [rand.(Normal(0, σ₀), size(x, 1)) rand.(Normal(), size(x, 1))]
    for i = 1:sum(Y .== 2)
        ind = ((i-1)*sims+1):(i*sims)
        if sum(I2(x[ind, :], z[ind, :], ε[ind, :]))!=0
            ll += log(1/sum(I2(x[ind, :], z[ind, :], ε[ind, :])))
            #ll += log(sum(I2(x[ind, :], z[ind, :], ε[ind, :])) / sims)
        else
            ll+=log(1/sims) #Act as though it happened at least once to avoid -infinity
        end
    end
    # Calculate log-likelihood for Y = 3 observations
    x, z = repeat(X[Y.==3, :], inner = [sims, 1]), repeat(Z[Y.==3, :], inner = [sims, 1])
    ε = [rand.(Normal(0, σ₀), size(x, 1)) rand.(Normal(), size(x, 1)) rand.(Normal(), size(x, 1))]
    for i = 1:sum(Y .== 3)
        ind = ((i-1)*sims+1):(i*sims)
        if sum(I3(x[ind, :], z[ind, :], ε[ind, :]))!=0
            #ll += log(1/sum(I3(x[ind, :], z[ind, :], ε[ind, :])))
            ll += log(sum(I3(x[ind, :], z[ind, :], ε[ind, :])) / sims)
        else
            ll+=log(1/sims) #Act as though it happened at least once to avoid -infinity
        end
    end
    # Calculate log-likelihood for Y = 4 observations
    x, z = repeat(X[Y.==4, :], inner = [sims, 1]), repeat(Z[Y.==4, :], inner = [sims, 1])
    ε = [rand.(Normal(0, σ₀), size(x, 1)) rand.(Normal(), size(x, 1)) rand.(Normal(), size(x, 1))]
    for i = 1:sum(Y .== 4)
        ind = ((i-1)*sims+1):(i*sims)
        if sum(I4(x[ind, :], z[ind, :], ε[ind, :]))!=0
            #ll += log(1/sum(I4(x[ind, :], z[ind, :], ε[ind, :])))
            ll += log(sum(I4(x[ind, :], z[ind, :], ε[ind, :])) / sims)
        else
            ll+=log(1/sims) #Act like it happened at least once to avoid -infinity
        end
    end
    return ll
end # accept-reject log-likelihood function
