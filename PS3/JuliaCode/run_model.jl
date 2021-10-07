using Distributed, SharedArrays

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

theme(:juno)
plot(out_res.val_fun[:,:, end])
plot!(out_res.val_fun[:,:, end-1])
plot!(out_res.val_fun[:,:, end-2])



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

plot(out_prim.a_grid, savings[:,1,20])
plot!(out_prim.a_grid, savings[:,2,20])

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
@time prim_noRisk, res_noRisk           = MarketClearing(use_Fortran=false, tol = 1e-3, i_risk = false);
@time prim_noRisk_noSS, res_noRisk_noSS = MarketClearing(use_Fortran=false, tol = 1e-3, ss = false, i_risk = false);
@time prim_exLab, res_exLab             = MarketClearing(use_Fortran=false, tol = 1e-3, exog_l = true);
@time prim_exLab_noSS, res_exLab_noSS   = MarketClearing(use_Fortran=false, tol = 1e-3, ss = false, exog_l = true);

# write results to table 1
open("PS3/Tables/table1.tex", "w") do io 
    write(io, string("\\begin{tabular}{|l|c|c|c|c|c|c|}\\hline",
    "&\\multicolumn{2}{c}{Benchmark Model} &\\multicolumn{2}{c}{No risk, \$z^L=z^H=0.5\$}",
    "&\\multicolumn{2}{c}{Exogenous labor, \$\\gamma=1\$}\\\\\\hline",
    "capital, \$K\$ & ", round(out_res.K, digits = 3), " & ", round(res_noSS.K, digits = 3), " & ", round(res_noRisk.K, digits = 3), " & ", round(res_noRisk_noSS.K, digits = 3), " & ",
    round(res_exLab.K, digits = 3), " & ", round(res_exLab_noSS.K, digits = 3), " \\\\\\hline",
    "labor, \$L\$ & ", round(out_res.L, digits = 3), " & ", round(res_noSS.L, digits = 3), " & ", round(res_noRisk.L, digits = 3), " & ", round(res_noRisk_noSS.L, digits = 3), " & ",
    round(res_exLab.L, digits = 3), " & ", round(res_exLab_noSS.L, digits = 3), " \\\\\\hline",
    "wage, \$w\$ & ", round(out_res.w, digits = 3), " & ", round(res_noSS.w, digits = 3), " & ", round(res_noRisk.w, digits = 3), " & ", round(res_noRisk_noSS.w, digits = 3), " & ",
    round(res_exLab.w, digits = 3), " & ", round(res_exLab_noSS.w, digits = 3), " \\\\\\hline",
    "interest, \$r\$ & ", round(out_res.r, digits = 3), " & ", round(res_noSS.r, digits = 3), " & ", round(res_noRisk.r, digits = 3), " & ", round(res_noRisk_noSS.r, digits = 3), " & ",
    round(res_exLab.r, digits = 3), " & ", round(res_exLab_noSS.r, digits = 3), " \\\\\\hline",
    "pension benefit, \$b\$ & ", round(out_res.b, digits = 3), " & ", round(res_noSS.b, digits = 3), " & ", round(res_noRisk.b, digits = 3), " & ", round(res_noRisk_noSS.b, digits = 3), " & ",
    round(res_exLab.b, digits = 3), " & ", round(res_exLab_noSS.b, digits = 3), " \\\\\\hline",
    "total welfare, \$W\$ & ", round(sum(out_res.val_fun.*out_res.F), digits = 3), " & ", round(sum(res_noSS.val_fun.*res_noSS.F), digits = 3), " & ", round(sum(res_noRisk.val_fun*res_noRisk.F), digits = 3), " & ",
    round(sum(res_noRisk_noSS.val_fun.*res_noRisk_noSS.F), digits = 3), " & ",round(sum(res_exLab.val_fun.*res_exLab.F), digits = 3), " & ", round(sum(res_exLab_noSS.val_fun.*res_exLab_noSS.F), digits = 3), " \\\\\\hline",
    "cv(wealth) & ", round(Lambda(put_prim, out_res, W), digits = 3), " & ", round(Lambda(prim_noSS, res_noSS, W), digits = 3), " & ", round(Lambda(prim_noRisk, res_noRisk, W), digits = 3),
    " & ", round(Lambda(prim_noRisk_noSS, res_noRisk_noSS, W), digits = 3), " & ", round(Lambda(prim_exLab, res_exLab, W), digits = 3), " & ", 
    round(Lambda(prim_exLab_noSS, res_exLab_noSS, W), digits = 3), " \\\\\\hline \\end{tabular}"))
end;
