# ------------------------------------------------------------------------------
# Author: Philip Coyle
# neogrowth_model_parallel.jl
# ------------------------------------------------------------------------------


## Structures
@with_kw struct Params
    β::Float64 = 0.99
    δ::Float64 = 0.025
    θ::Float64 = 0.36

    tol::Float64 = 1e-4
    maxit::Int64 = 10000
end

@with_kw struct Shocks
    Zg::Float64 = 1.25
    Zb::Float64 = 0.2

    # Transition Probabilities
    # Pyx = Pr(Z' = Zy | Z = Zx) x∈{g,b} & y∈{g,b}
    Pgg::Float64 = 0.977
    Pbg::Float64 = 1 - Pgg
    Pbb::Float64 = 0.926
    Pgb::Float64 = 1 - Pbb
    Tmat::Array{Float64,2} = [
        Pgg Pbg
        Pgb Pbb
    ]
end

@with_kw struct Grids
    Zg::Float64 = 1.25
    Zb::Float64 = 0.2

    k_lb::Float64 = 0.01
    k_ub::Float64 = 45.0
    n_k::Int64 = 1000
    k_grid::Array{Float64,1} = range(k_lb, stop = k_ub, length = n_k)

    n_Z::Int64 = 2
    Z_grid::Array{Float64,1} = [Zg, Zb]
end

mutable struct PolFuncs
    pf_c::SharedArray{Float64,2}
    pf_k::SharedArray{Float64,2}
    pf_v::SharedArray{Float64,2}
end


## Functions
function solve_model()
    P = Params()
    S = Shocks()
    G = Grids()

    @unpack n_k, n_Z = G
    # Initial Guess
    pf_c = SharedArray{Float64}(zeros(n_k, n_Z))
    pf_k = SharedArray{Float64}(zeros(n_k, n_Z))
    pf_v = SharedArray{Float64}(zeros(n_k, n_Z))
    PFs = PolFuncs(pf_c, pf_k, pf_v)

    converged = 0
    it = 1
    while converged == 0 && it < P.maxit
        @unpack pf_v, pf_k, pf_c = PFs

        pf_c_up, pf_k_up, pf_v_up = Bellman(P, S, G, PFs)

        diff_v = maximum(abs.(pf_v_up - pf_v))
        diff_k = maximum(abs.(pf_k_up - pf_k))
        diff_c = maximum(abs.(pf_c_up - pf_c))

        max_diff = diff_v + diff_k + diff_c

        if mod(it, 50) == 0 || max_diff < P.tol
            println(" ")
            println("*************************************************")
            println("AT ITERATION = ", it)
            println("MAX DIFFERENCE = ", max_diff)
            println("*************************************************")

            if max_diff < P.tol
                converged = 1
            end
        end
        # Update the policy functions
        PFs = PolFuncs(pf_c_up, pf_k_up, pf_v_up)
        it = it + 1
    end

    return G, PFs
end

function Bellman(P::Params, S::Shocks, G::Grids, PFs::PolFuncs)
    @unpack β, δ, θ = P
    @unpack Zg, Zb, Pgg, Pbg, Pbb, Pgb, Tmat = S
    @unpack n_k, k_grid, n_Z, Z_grid = G
    @unpack pf_c, pf_k, pf_v = PFs

    # To make updating work
    pf_k_up = SharedArray{Float64}(zeros(n_k, n_Z))
    pf_c_up = SharedArray{Float64}(zeros(n_k, n_Z))
    pf_v_up = SharedArray{Float64}(zeros(n_k, n_Z))

    # @sync @distributed for (i_Z, Z) in enumerate(Z_grid)
    @sync @distributed for i_Z = 1:n_Z
                           Z = Z_grid[i_Z]
                           Pr = Tmat[i_Z, :]
                           for (i_k, k_today) in enumerate(k_grid)
                               # Must be defined outside loop.
                               v_today = log(0)
                               c_today = log(0)
                               k_tomorrow = log(0)

                               y_today = Z * k_today^θ

                               # Find optimal investment/consumption given capital level today
                               for (i_kpr, k_temp) in enumerate(k_grid)
                                   c_temp = y_today + (1 - δ) * k_today - k_temp
                                   v_tomorrow = Pr[1] * pf_v[i_kpr, 1] + Pr[2] * pf_v[i_kpr, 2]
                                   if c_temp < 0
                                       v_temp = log(0) + β * v_tomorrow
                                   else
                                       v_temp = log(c_temp) + β * v_tomorrow
                                   end

                                   if v_temp > v_today
                                       v_today = v_temp
                                       c_today = c_temp
                                       k_tomorrow = k_temp
                                   end
                               end

                               # Update PFs
                               pf_k_up[i_k, i_Z] = k_tomorrow
                               pf_c_up[i_k, i_Z] = c_today
                               pf_v_up[i_k, i_Z] = v_today
                           end
                       end

    return pf_c_up, pf_k_up, pf_v_up
end

function plot_pfs(dir::String, G::Grids, PFs::PolFuncs)
    @unpack k_grid = G
    @unpack pf_c, pf_k, pf_v = PFs

    pf1 = plot(k_grid,pf_v[:,1],legend = false,color=:black, label = "Good State",lw = 2);
    pf_1 = plot!(pf1,k_grid,pf_v[:,2],title="Value Function",legend = true,color=:blue, label = "Bad State",lw = 2);

    pf2 = plot(k_grid,pf_k[:,1],legend = false,color=:black, lw = 2);
    pf_2 = plot!(pf2,k_grid,pf_k[:,2],title="Capital Investment",legend = false,color=:blue, lw = 2);

    # pf3 = plot(k_grid,pf_c[:,1],legend = false,color=:black, lw = 2);
    # pf_3 = plot!(pf3,k_grid,pf_c[:,2],title="Consumption",legend = false,color=:blue, lw = 2);

    pf3 = plot(k_grid,pf_k[:,1] - k_grid,legend = false,color=:black, lw = 2);
    pf_3 = plot!(pf3,k_grid,pf_k[:,2] - k_grid, title="Net Capital Inv.",legend = false,color=:blue, lw = 2);

    pf = plot(pf_1,pf_2,pf_3,layout=(1,3),size = (600,400)) #Size can be adjusted so don't need to mess around with 'blank space'
    xlabel!("Initial Capital Stock")
    # savefig(dir * "Q7_PFs.pdf")
end
