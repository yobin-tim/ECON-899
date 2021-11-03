#==
    This file defines functions used in JF's PS1
==#

# Calculate log-likelihood at β
function likelihood(β, Y, X)

    sum(Y.*log.(exp.(X*β) ./ (1 .+ exp.(X*β))) +
        (1 .- Y).*log.(1 ./ (1 .+ exp.(X*β))))

end # log-likelihood function

# calculate the log-likelihood score, given β
function score(β, Y, X)

    sum((Y .- (exp.(X*β_0) ./ (1 .+ exp.(X*β_0)))) .* X, dims = 1)

end # end log-likelihood score

# Calculate the Hessian matrix given β
function Hessian(X, β)

    A = (exp.(X*β) ./ (1 .+ exp.(X*β))) .*
        (1 ./ (1 .+ exp.(X*β_0)))

    B = zeros(size(X,2), size(X,2), size(X,1))

    for i = 1:size(X,1)

        B[:,:,i] = A[i] .* X[i,:] * transpose(X[i,:])

    end

    dropdims(sum(B, dims = 3), dims = 3)

end # Hessian matrix

# Define the Newton convergence algorithm
function Newton(Y, X; β₀ = [-1.0; ones(size(X, 2), 1)], err::Float64 = 100, tol::Float64 = 10e-8)

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