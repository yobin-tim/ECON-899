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
        H = H .+ (A[i] .* X[i,:] * transpose(X[i,:]))
    end

    return -H
end # Hessian matrix

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
function NewtonAlg(Y, X; β₀::Matrix{Float64} = [-1.0; ones(size(X, 2), 1)], err::Float64 = 100.0, tol::Float64 = 1e-8)
    β_out=0
    while err > tol

        # update β
        β_out = β₀ - inv(Hessian(X, β₀))*transpose(score(β₀, Y, X))

        # calculate error and update β₀
        err = maximum(abs.(β_out - β₀))
        β₀ = copy(β_out)

    end # err > tol loop

    # return converged β
    return β_out
end # Newton's algorithm
