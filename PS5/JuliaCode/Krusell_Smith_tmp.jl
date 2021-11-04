using Distributed, SharedArrays, Interpolations, Optim

workers()
addprocs(1)

@Distributed.everywhere include("./HelpfulFunctions.jl")

prim, grid, shock, res = Initialize();

ℇ_idio, ℇ_agg = draw_shocks(shock, prim.N, prim.T)


####################################################
####################################################

## Return pf_v and pf_k in res.
VFI(prim, grid, shock, res)

## time = 1 to get \bar_K_2
grid.K_grid[4] ## 11.5 ~ K_SS
z = Int.(ℇ_agg)[1] ## Aggrigate state time = 1 
ϵ = Int.(ℇ_idio)[:,1] ## Idio state time = 1
k1 = zeros(5000)

for i = 1:5000
    k1[i] = res.pf_k[4,ϵ[i],4, z] ## Calculate each persons saving.
end

K_2 = mean(k1) ## Aggregate Capital = 3.59

## There is no values in K grid (3.59)
## in k grid, there is 3.00085 [4] and 4.0008 [5].
## How to define K_3?

z = Int.(ℇ_agg)[2]
ϵ = Int.(ℇ_idio)[:,2]
k2 = zeros(5000)

for i = 1:5000
    k2[i] = res.pf_k[??,ϵ[i],??, z] ## How to define the index?
end

K_3 = mean(k2)
#######################################################


