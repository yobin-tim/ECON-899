#==
    This file defines functions used in JF's PS1
==#
using Optim

# Calculate log-likelihood at β
function likelihood(β, Y, X)

    X = [ones(size(X, 1), 1) X] # add constant to X

    return sum(Y.*log.(exp.(X*β) ./ (1 .+ exp.(X*β))) +
        (1 .- Y).*log.(1 ./ (1 .+ exp.(X*β))))
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
        H = H .+ A[i] .* X[i,:] * transpose(X[i,:])
    end

    return H
end # Hessian matrix

#Calculate FirstDerivate
function ∂F(β,Y,X;h=1e-7)
    ∂=zeros(length(β))
    for ii=1:length(β)
        hi=zeros(length(β))
        hi[ii]=copy(h)
        ∂[ii]=(likelihood(β.+h,Y,X)-likelihood(β,Y,X))/h
    end
    return transpose(∂)
end
# Define the Newton convergence algorithm
function NewtonAlg(Y, X; β₀::Matrix{Float64} = [-1.0; ones(size(X, 2), 1)], err::Float64 = 100, tol::Float64 = 10e-8)

    while err > tol

        # update β
        β = β₀ - inv(Hessian(X, β₀))*score(β₀, Y, X)

        # calculate error and update β₀
        err = maximum(abs.(β - β₀))
        β₀ = β

    end # err > tol loop

    # return converged β
    return β

end # Newton's algorithm
