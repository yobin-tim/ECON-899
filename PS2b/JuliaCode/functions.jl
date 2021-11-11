#==
    This file defines functions used in JF's PS2
==#
using Optim, Distributions

# structure of model parameters
mutable struct Parameters
    α₀::Float64
    α₁::Float64
    α₂::Float6
    β::Array{Float64}
    γ::Array{Float64}
    ρ::Float64
end # parameters struct

# Calculate log-likelihood at β
function likelihood(Y, X, Z, W1, W2, param::Parameters)
    # separate weights and nodes from W1 and W2
    u = W1[:, 1]; w = W1[:, 2]
    μ0 = W2[:, 1]; μ1 = W2[:, 2]; ω = W2[:, 3]

    # unpack model parametersw
    @unpack α₀, α₁, α₂, β, γ, ρ = param

    # Calculate σ₀ and σ₀²
    σ₀² = 1/(1-ρ)^2
    σ₀  = 1/(1-ρ)

    X = [ones(size(X, 1), 1) X] # add constant to X

    # calculate likelihood for each level of Y 
    Y1 = Y[Y == 1]
    Y2 = Y[Y == 2]
    Y3 = Y[Y == 3]
    Y4 = Y[Y == 4]

    # map integral bounds to [0, 1]
    b₀ = (x) -> log(x) + α₀ - X*β - Z*γ
    b₁ = (x) -> log(x) + α₁ - X*β - Z*γ

    # Y = 1 likelihood:
    L1 = cdf(Normal(), (-α₀ - X*β - Z*γ)/σ₀²)
    L2 = sum(w.*(cdf.(Normal(), (-ρ).*b₀.(u) .- (α₁ + X*β + Z*γ))./σ₀).*pdf.(Normal(), b₀.(u./σ₀)))
    L3 = sum(ω.*((cdf.(Normal(), (-ρ).*b₁(μ₁) .- (α₂ + X*β + Z*γ)))./σ₀).*pdf.(Normal(), b₀.(μ₀./σ₀).*pdf.(Normal(), b₁(μ₁).-ρ.*b₀(μ₀))))
    L4 = sum(ω.*(cdf.(Normal(), (-ρ).*b₁(μ₁) .+ (α₂ + X*β + Z*γ))./σ₀).*pdf.(Normal(), b₁(μ₁./σ₀).*pdf.(Normal(), b₁(μ₁).-ρ.*b₀(μ₀))))

    return sum((Y .== 1).*log.(L1) + (Y .== 2).*log.(L2) + (Y .== 3).*log.(L3) + (Y .== 4).*log.(L4))
end # log-likelihood function

# calculate the log-likelihood score, given β
function score(β, Y, X)

    X = [ones(size(X, 1), 1) X] # add constant to X

    return sum((Y .- (exp.(X*β) ./ (1 .+ exp.(X*β)))) .* X, dims = 1)

end # end log-likelihood score

# Calculate the Hessian matrix given β
function Hessian(X, β)

    X = [ones(size(X, 1), 1) X] # add constant to X

    A = (exp.(X*β) ./ (1 .+ exp.(X*β))) .*
        (1 ./ (1 .+ exp.(X*β)))

    #==
    B = zeros(size(X,2), size(X,2), size(X,1))

    for i = 1:size(X,1)

        B[:,:,i] = A[i] .* X[i,:] * transpose(X[i,:])

    end

    dropdims(sum(B, dims = 3), dims = 3)
    ==#

    # Alternative method (saves memory):
    H = 0;
    for i = 1:size(X,1)
        H = H .+ (A[i] .* X[i,:] * transpose(X[i,:]))
    end

    return -H
end # Hessian matrix

####################################################################
#Calculate First Derivate numerically
function ∂F(β,Y,X;h=1e-5)
    ∂=zeros(length(β))
    for ii=1:length(β)
        hi=zeros(length(β))
        hi[ii]=copy(h)
        ∂[ii]=(likelihood(β.+h,Y,X)-likelihood(β,Y,X))/h
    end
    return transpose(∂)
end

function score_num(β,Y,X;h=1e-5)

    partial = zeros(length(β))
    for i =1:length(β)
        β1=copy(β)
        β1[i] += h 
        partial[i]=(likelihood(β1,Y,X)-likelihood(β,Y,X))/h
    end
    return transpose(partial)
end

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
