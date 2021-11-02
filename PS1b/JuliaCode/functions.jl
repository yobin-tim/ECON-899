using CSV
using DataFrames

df = DataFrame(CSV.File("../data/mortgage.csv"))

#df[!, :i_25] = df[!, :i_open_year2] .- df[!, :i_open_year5]
#I think i_open_year2-i_open_year5 is Stata notation for
#"include variables i_open_year2 through i_open_year5"
X = df[!,[:i_large_loan,:i_medium_loan,:rate_spread,
          :i_refinance,:age_r,:cltv,:dti, :cu,
          :first_mort_r,:score_0,:score_1, :i_FHA,
          :i_open_year2,:i_open_year3, :i_open_year4,
          :i_open_year5]] |> Matrix

Y = df[!, :i_close_first_year] #|> Matrix


function likelihood(β, Y, X)

    sum(Y.*log.(exp.(X*β) ./ (1 .+ exp.(X*β))) +
        (1 .- Y).*log.(1 ./ (1 .+ exp.(X*β))))

end

function score(β, Y, X)

    sum((Y .- (exp.(X*β_0) ./ (1 .+ exp.(X*β_0)))) .* X, dims = 1)

end

function Hessian(X, β)

    A = (exp.(X*β) ./ (1 .+ exp.(X*β))) .*
        (1 ./ (1 .+ exp.(X*β_0)))

    B = zeros(size(X,2), size(X,2), size(X,1))

    for i = 1:size(X,1)

        B[:,:,i] = A[i] .* X[i,:] * transpose(X[i,:])

    end

    dropdims(sum(B, dims = 3), dims = 3)

end
