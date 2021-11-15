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

    # map integral bounds to [0, 1]
    b₀ = (a, x, z) -> log.(a) .+ (α₀ + dot(x, β) + dot(z, γ))
    b₁ = (a, x, z) -> log.(a) .+ (α₁ + dot(x, β) + dot(z, γ))

    # save indexes for each outcome
    indvec = 1:size(Y, 1)
    ind₁ = indvec[Y .== 1]; ind₂ = indvec[Y .== 2];
    ind₃ = indvec[Y .== 3]; ind₄ = indvec[Y .== 4];

    # per-observation likelihood:
    L1 = (x, z) -> log(cdf(Normal(), (-α₀ - dot(x, β) - dot(z, γ))/σ₀)) # I think this should just be σ₀ not σ₀²
    L2 = (x, z) -> log(sum(w.*(cdf.(Normal(), (-ρ)*b₀(u, x, z) .- (α₁ .+ dot(x, β) +
        dot(z, γ)))./σ₀).*pdf.(Normal(), b₀(u./σ₀, x, z)) ./ u))
    L3 = (x, z) -> log(sum(ω.*((cdf.(Normal(), (-ρ)*b₁(μ₁, x, z) .- (α₂ .+ dot(x, β) +
        dot(z, γ))))./σ₀).*pdf.(Normal(), b₀(μ₀./σ₀, x, z)).*pdf.(Normal(), b₁(μ₁, x, z).-
        ρ*b₀(μ₀, x, z)) ./ (μ₀ .* μ₁)))
    L4 = (x, z) -> log(sum(ω.*(cdf.(Normal(), (ρ)*b₁(μ₁, x, z) .+ (α₂ + dot(x, β) +
        dot(z, γ)))./σ₀).*pdf.(Normal(), b₁(μ₁./σ₀, x, z)).*pdf.(Normal(), b₁(μ₁, x, z) .-
        ρ*b₀(μ₀, x, z)) ./ (μ₀ .* μ₁)))

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
function GHKLL(Y, X, Z, θ; sims = 100, k = maximum(Y))

    # unpack model parameters
    param = ModelParameters(θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
    @unpack α₀, α₁, α₂, β, γ, ρ = param
    σ₀  = 1/(1-ρ)

    #==
    TODO: figure out how to do this; I'm lost -Danny
    # Draw ν₁
    Φ₁ = truncated.(Normal(), -Inf, -X*β)
    η₁ = rand.(Φ₁, sims)

    # Draw ν₂
    Φ₂ = truncated.(Normal(), -Inf, -X*β)
    ==#
    #=
        Maybe it is something along the lines of what I have below? I will try to do more after lunch
        Do we do separate random draws for each i? That is what I will do here ~Ryan
    =#
    ll=0
    for i=1:size(Y,1)
        ϵ_draws=zeros(sims,minimum(Y[i],3))
        if Y[i]>1 #Need to draw from a distribution which won't make the borrower repay in period 1
            ϵ_draws[:,1]=rand(truncated(Normal(0,σ₀),-Inf,-α₀-X[i,:]*β-Z[i,:]*γ),sims)
        elseif Y[i]==1 #Need to draw from a distribution which will make the borrower repay in period 1
            ϵ_draws[:,1]=rand(truncated(Normal(0,σ₀),-α₀-X[i,:]*β-Z[i,:]*γ,Inf),sims)
        end
        if Y[i]>2 #Need to draw from a distribution which won't make the borrower repay in period 2
            ϵ_draws[:,2]=rand(truncated(Normal(0,σ₀),-Inf,-α₀-X[i,:]*β-Z[i,:]*γ-ρ*ϵ_draws[:,1]),sims)
        elseif Y[i]==2 #Need to draw from a distribution which will make the borrower repay in period 2
            ϵ_draws[:,2]=rand(truncated(Normal(0,σ₀),-α₀-X[i,:]*β-Z[i,:]*γ-ρ*ϵ_draws[:,1],Inf),sims)
        end
        if Y[i]>3 #Need to draw from a distribution which won't make the borrower repay in period 3
            ϵ_draws[:,3]=rand(truncated(Normal(0,σ₀),-Inf,-α₀-X[i,:]*β-Z[i,:]*γ-(ρ^(2))*ϵ_draws[:,1]-ρ*ϵ_draws[:,2]),sims)
        elseif Y[i]==3 #Need to draw from a distribution which will make the borrower repay in period 4
            ϵ_draws[:,3]=rand(truncated(Normal(0,σ₀),-α₀-X[i,:]*β-Z[i,:]*γ-(ρ^(2))*ϵ_draws[:,1]-ρ*ϵ_draws[:,2],Inf),sims)
        end

    end


end # quadrature log-likelihood function






######################################################################

#Calculate the Hessian numerically
function Find_H_num(β,Y,X;h=1e-5)
    H_num=zeros(length(β),length(β))
    d=1
    for i1=1:length(β)
        for i2=d:length(β)
            h1=zeros(length(β))
            h2=copy(h1)
            h1[i1], h2[i2] = copy(h),copy(h)
            #This formula was taken from http://www.holoborodko.com/pavel/2014/11/04/computing-mixed-derivatives-by-finite-differences/
            #H_num[i1,i2]=(likelihood(β.-h1.-h2,Y,X)+likelihood(β.+h1.+h2,Y,X)+
            #    likelihood(β.+h1.-h2,Y,X)+likelihood(β.-h1.+h2,Y,X))/(4*(h^2))
            #Alternate, more accurate formula also from the above link
            H_num[i1,i2]=( 8*(likelihood(β.+h1.-2 .*h2,Y,X)+likelihood(β.+2 .* h1.-h2,Y,X)+
                            likelihood(β.-2 .*h1.+h2,Y,X)+likelihood(β.-h1.+2 .*h2,Y,X))-
                            8*(likelihood(β.-h1.-2 .*h2,Y,X)+likelihood(β.-2 .* h1.-h2,Y,X)+
                            likelihood(β.+h1.+ 2 .*h2,Y,X)+likelihood(β.+ 2 .* h1.+h2,Y,X))-
                            (likelihood(β.+2 .*h1.-2 .*h2,Y,X)+likelihood(β.-2 .* h1.+ 2 .*h2,Y,X)-
                            likelihood(β.-2 .*h1.- 2 .*h2,Y,X)-likelihood(β.+ 2 .*h1.+2 .*h2,Y,X))+
                            64*(likelihood(β.-h1.-h2,Y,X)+likelihood(β.+h1.+h2,Y,X)-
                            likelihood(β.+h1.-h2,Y,X)-likelihood(β.-h1.+h2,Y,X)))/(144*(h^2))
        end
        d+=1
    end
    #Exploit Hessian symmetry to find the remaining entries
    d=length(β)-1
    for i1=length(β):-1:2
        for i2=d:-1:1
            H_num[i1,i2]=H_num[i2,i1]
        end
        d-=1
    end
    return H_num
end
# Define the Newton convergence algorithm
function NewtonAlg(Y, X; β₀::Matrix{Float64} = [-1.0; ones(size(X, 2), 1)],
        err::Float64 = 100.0, tol::Float64 = 1e-32, sk::Float64=1e-7)
    β_out=0
    iter=1;
    print("="^35,"\n","Newton's Method","\n")
    while err > tol

        # update β
        β_out = β₀ - sk*inv(Hessian(X, β₀))*transpose(score(β₀, Y, X))
        #If you have made β_out NaN, you've gone too far
            while isnan((transpose(β_out)ones(size(X, 2)+1, 1))[1])
                sk=sk/10;
                β_out = β₀ - sk*inv(Hessian(X, β₀))*transpose(score(β₀, Y, X))
            end
        # calculate error and update β₀
        err_new = maximum(abs.(β_out - β₀))
        β₀ = copy(β_out)
        if iter % 5==0
            println("Newton Iteration $(iter) with error $(err_new)")
        end
        iter+=1
        #Update sk depending on whether things are going well or not
            if err_new<err
                sk=sk*2;
            else
                sk=sk/10
            end
            err=copy(err_new)
    end # err > tol loop

    # return converged β
    print("\nThe Newton algorithm takes")
    return β_out
end # Newton's algorithm
