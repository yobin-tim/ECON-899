@time begin
using Distributed, SharedArrays, NaNMath, JLD

#add processes
workers()
addprocs(2)


#@Distributed.everywhere include("./PS3/JuliaCode/conesa_kueger.jl");
@Distributed.everywhere include("./conesa_kueger.jl");

#prim, res = Initialize(); #=
#@time V_ret(prim, res);
#@time V_workers(prim, res);
#@time V_Fortran(res.r, res.w, res.b);
#@time SteadyStateDist(prim, res); =#

#agridf, consumption = V_Fortran(res.r, res.w, res.b);
#=
agridf
prim.a_grid
hcat(prim.a_grid, agridf, prim.a_grid - agridf)
=#
@time out_prim, out_res = MarketClearing(use_Fortran=false, tol = 1e-3);

using Plots, LaTeXStrings

theme(:vibrant)
default(fontfamily = "Computer Modern", framestyle=:box) #LaTex Style
plot(out_res.val_fun[:,:, end])
plot!(out_res.val_fun[:,:, end-1])
plot!(out_res.val_fun[:,:, end-2])

# Exercise 1
plot(out_prim.a_grid, out_res.val_fun[:,:, 50],
     title = "Value Function at age 50",
     xlabel = "Assets",
     ylabel = "Value Function",
     color=:black,
     legend = false,
     lw = 2)
savefig("../Figures/value_function50.pdf")


plot(out_prim.a_grid, out_res.val_fun[:, 1, end])
plot!(out_prim.a_grid, out_res.val_fun[:, 2, end])
plot(out_prim.a_grid, out_res.val_fun[:, 1, 20])
plot!(out_prim.a_grid, out_res.val_fun[:, 2,20])

plot(out_prim.a_grid, out_res.pol_fun[:,1,20])
plot!(out_prim.a_grid,out_res.pol_fun[:,2,20])

plot(out_prim.a_grid, out_res.val_fun[:, 1, 34])
plot!(out_prim.a_grid, out_res.val_fun[:, 2, 34])

# Calculate savings
savings = out_res.pol_fun[:,:,:] .- out_prim.a_grid;

plot(out_prim.a_grid, savings[:,1,20],
     title = "Saving Function at age 20",
     xlabel = "Assets",
     ylabel = "Saving Functions",
     label = "High productivity",
     color=:black,
     lw = 2,
     legend=:topright)
plot!(out_prim.a_grid, savings[:,2,20],
      label = "Low productivity",
      color=:black,
      line=:dash)
savefig("../Figures/savings_20.pdf")


#Plotting Policy Functions
    #Find a_hat (The asset point beyond which everyone dissaves)
    function Find_a_hat()
        a_hat=zeros(out_prim.nZ, out_prim.N_final)
        for zi=1:out_prim.nZ, ni=1:out_prim.N_final
            for ai=1:out_prim.nA
                if out_res.pol_fun[ai,zi,ni]<=out_prim.a_grid[ai]
                    a_hat[zi,ni]=out_prim.a_grid[ai]
                    break
                end
            end
        end
        return a_hat
    end
    a_hat=Find_a_hat()
    plot(out_prim.a_grid, out_res.pol_fun[:,1,20])
    plot!(out_prim.a_grid, out_res.pol_fun[:,2,20])
    vline!([a_hat[1,20]], label=L"\hat{a}")

    plot(1:out_prim.N_final,a_hat[1,:],
        xlabel="Age", ylabel=L"\hat{a}",label="High Shock")
    plot!(1:out_prim.N_final,a_hat[2,:], title="Maximum Assets at which an Agent Saves
over the Life Cycle",
        xlabel="Age", ylabel=L"\hat{a}",label="Low Shock")

plot!(out_prim.a_grid, out_prim.a_grid)

plot(savings[500,1,:])
plot!(savings[500,2,:])
plot!

# graph distributions

# accross assets for workers and retirees
a_dist = sum(out_res.F, dims = 2:3)
plot(out_prim.a_grid, a_dist[:, 1, 1])

plot(out_prim.a_grid,out_res.F[:,1,50])

# conduct policy experiments
@time prim_noSS, res_noSS               = MarketClearing(use_Fortran=false, tol = 1e-3, ss = false);

    ## for PS4 
    save("../../PS4/Data/Initial_Conditions.jld",
         "Γ_0", out_res.F,
         "V_0", out_res.val_fun,
         "K_θ", out_res.K,
         "K", res_noSS.K,
         "L_θ", out_res.L,
         "L", res_noSS.L,
         "V", res_noSS.val_fun,
         "a", res_noSS.pol_fun,
         "l", res_noSS.l_fun
         )

@time prim_noRisk, res_noRisk           = MarketClearing(use_Fortran=false, tol = 1e-2, i_risk = false);
@time prim_noRisk_noSS, res_noRisk_noSS = MarketClearing(use_Fortran=false, tol = 1e-2, ss = false, i_risk = false);
@time prim_exLab, res_exLab             = MarketClearing(use_Fortran=false, tol = 1e-3, exog_l = true);
@time prim_exLab_noSS, res_exLab_noSS   = MarketClearing(use_Fortran=false, tol = 1e-3, ss = false, exog_l = true);

# write results to table 1
open("../Tables/table1.tex", "w") do io 
    write(io, string(
        "\\begin{table}", 
        "\n \\centering",
        "\n \\caption{\\label{tab1} Results of policy experiments}",
        "\\begin{tabular}{lcccccc}",
        "\n \\hline",
        "\n \\hline",
        "\n &\\multicolumn{2}{c}{Benchmark Model} &\\multicolumn{2}{c}{No risk, \$z^L=z^H=0.5\$}&\\multicolumn{2}{c}{Exogenous labor, \$\\gamma=1\$}\\\\",
        "\n \\hline",         
        "\n capital, \$K\$ & ", round(out_res.K, digits = 3), " & ", round(res_noSS.K, digits = 3), " & ", round(res_noRisk.K, digits = 3), " & ", round(res_noRisk_noSS.K, digits = 3), " & ", round(res_exLab.K, digits = 3), " & ", round(res_exLab_noSS.K, digits = 3), " \\\\",
        "\n labor, \$L\$ & ", round(out_res.L, digits = 3), " & ", round(res_noSS.L, digits = 3), " & ", round(res_noRisk.L, digits = 3), " & ", round(res_noRisk_noSS.L, digits = 3), " & ", round(res_exLab.L, digits = 3), " & ", round(res_exLab_noSS.L, digits = 3), " \\\\",
        "\n wage, \$w\$ & ", round(out_res.w, digits = 3), " & ", round(res_noSS.w, digits = 3), " & ", round(res_noRisk.w, digits = 3), " & ", round(res_noRisk_noSS.w, digits = 3), " & ", round(res_exLab.w, digits = 3), " & ", round(res_exLab_noSS.w, digits = 3), " \\\\",
        "\n interest, \$r\$ & ", round(out_res.r, digits = 3), " & ", round(res_noSS.r, digits = 3), " & ", round(res_noRisk.r, digits = 3), " & ", round(res_noRisk_noSS.r, digits = 3), " & ", round(res_exLab.r, digits = 3), " & ", round(res_exLab_noSS.r, digits = 3), " \\\\",
        "\n pension benefit, \$b\$ & ", round(out_res.b, digits = 3), " & ", round(res_noSS.b, digits = 3), " & ", round(res_noRisk.b, digits = 3), " & ", round(res_noRisk_noSS.b, digits = 3), " & ", round(res_exLab.b, digits = 3), " & ", round(res_exLab_noSS.b, digits = 3), " \\\\",
        "\n total welfare, \$W\$ & ", round(NaNMath.sum(out_res.val_fun.*out_res.F), digits = 3), " & ", round(NaNMath.sum(res_noSS.val_fun.*res_noSS.F), digits = 3), " & ", round(NaNMath.sum(res_noRisk.val_fun.*res_noRisk.F), digits = 3), " & ", round(NaNMath.sum(res_noRisk_noSS.val_fun.*res_noRisk_noSS.F), digits = 3), " & ",round(NaNMath.sum(res_exLab.val_fun.*res_exLab.F), digits = 3), " & ", round(NaNMath.sum(res_exLab_noSS.val_fun.*res_exLab_noSS.F), digits = 3), " \\\\",
        "\n cv(wealth) & \\textemdash & ", round(Lambda(prim_noSS, res_noSS, out_res.val_fun), digits = 3), " & \\textemdash & ", round(Lambda(prim_noRisk_noSS, res_noRisk_noSS, res_noRisk.val_fun), digits = 3), " & \\textemdash & ", round(Lambda(prim_exLab_noSS, res_exLab_noSS, res_exLab.val_fun), digits = 3), " \\\\",
        "\n \\hline",
        "\n \\end{tabular}",
        "\n \\end{table}"))
end;


# Exercise 3
c = Lambda2(prim_noSS, res_noSS, out_res.val_fun);

plot(c[1], c[2],
     title = "",
     xlabel = "Age",
     ylabel = "λ",
     color=:black,
     legend = false,
     lw = 2)
savefig("../Figures/lambda.pdf")
end
