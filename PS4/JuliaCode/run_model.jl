using Distributed, SharedArrays, JLD

#add processes
workers()
addprocs(2)


@Distributed.everywhere include("./PS4/JuliaCode/conesa_kueger.jl");
#@Distributed.everywhere include("./conesa_kueger.jl");

#Exercise 1: (Problem Set 4)
@time out_K_path,out_Ft= TransitionPath(TrySaveMethod=false,Experiment=1)
#Exercise 2: (Problem Set 4)
@time out_K_path_Exp2,out_Ft_Exp2= TransitionPath(TrySaveMethod=false,Experiment=2)

using Plots, LaTeXStrings

theme(:juno)
plot(1:151,out_K_path[:], ylabel="Aggregate Capital",
    xlabel="Time",label="Excercise 1")
    plot!(1:151,out_K_path_Exp2[:], legend=:bottomright, ylabel="Aggregate Capital",
        xlabel="Time",label="Excercise 2")
    savefig("./PS4/Document/Figures/ComparingTransitions.png")

function Exercise1Prob2(kpath; α=.36,δ=.06,N_final=66,J_R=46, TransitionNo=151)
    #Recalculating Aggregate Labor in the Inelastic case
        L=0.7543 #This was the converged value of inelastic labor supply
    #Functions for interest rate and wages
        r_mkt   ::Function          = (K, L) -> α*(K^(α-1))*(L^(1-α)) - δ
        w_mkt   ::Function          = (K, L) -> (1-α)*(K^α)*(L^(-α))
    r_trans=[r_mkt(k,L) for k in kpath]
    w_trans=[w_mkt(k,L) for k in kpath]
    #Plotting Aggregate Capital
        plot(1:TransitionNo,out_K_path[:], ylabel="Aggregate Capital",
            xlabel="Time",legend=false)
        savefig("./PS4/Document/Figures/PathOfAggregateCapital.png")
        plot(1:TransitionNo,r_trans[:], ylabel="Interest Rate",
            xlabel="Time",legend=false)
        savefig("./PS4/Document/Figures/PathOfInterestRate.png")
        plot(1:TransitionNo,w_trans[:], ylabel="Wages",
            xlabel="Time",legend=false)
        savefig("./PS4/Document/Figures/PathOfWages.png")
end
Exercise1Prob2(out_K_path)
#= Old plots for problem set 3
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
        xlabel="Age", ylabel=L"\hat{a}",label="High Schock")
    plot!(1:out_prim.N_final,a_hat[2,:], title="Maximum Assets at which an Agent Saves
over the Life Cycle",
        xlabel="Age", ylabel=L"\hat{a}",label="Low Schock")

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
    write(io, string(L"\begin{tabular}{|l|c|c|c|c|c|c|}\hline",
    L"&\multicolumn{2}{c}{Benchmark Model} &\multicolumn{2}{c}{No risk, $z^L$=$z^H=0.5$}",
    L"&\multicolumn{2}{c}{Exogenous labor, $\gamma=1$}\\\hline",
    L"capital, $K$ & ", out_res.K, " & ", res_noSS.K, " & ", res_noRisk.K, " & ", res_noRisk_noSS.K, " & ",
    res_exLab.K, " & ", res_exLab_noSS.K, " \\\hline",
    L"labor, $L$ & ", out_res.L, " & ", res_noSS.L, " & ", res_noRisk.L, " & ", res_noRisk_noSS.L, " & ",
    res_exLab.L, " & ", res_exLab_noSS.L, " \\\hline",
    L"wage, $w$ & ", out_res.w, " & ", res_noSS.w, " & ", res_noRisk.w, " & ", res_noRisk_noSS.w, " & ",
    res_exLab.w, " & ", res_exLab_noSS.w, " \\\hline",
    L"interest, $r$ & ", out_res.r, " & ", res_noSS.r, " & ", res_noRisk.r, " & ", res_noRisk_noSS.r, " & ",
    res_exLab.r, " & ", res_exLab_noSS.r, " \\\hline",
    L"pension benefit, $b$ & ", out_res.b, " & ", res_noSS.b, " & ", res_noRisk.b, " & ", res_noRisk_noSS.b, " & ",
    res_exLab.b, " & ", res_exLab_noSS.b, " \\\hline",
    L"total welfare, $W$ & ", sum(out_res.val_fun.*out_res.F), " & ", res_noSS.b, " & ", res_noRisk.b, " & ", res_noRisk_noSS.b, " & ",
    res_exLab.b, " & ", res_exLab_noSS.b, " \\\hline",
    L"cv(wealth) & ", Lambda(put_prim, out_res, W), " & ", Lambda(prim_noSS, res_noSS, W), " & ", Lambda(prim_noRisk, res_noRisk, W),
    " & ", Lambda(prim_noRisk_noSS, res_noRisk_noSS, W), " & ", Lambda(prim_exLab, res_exLab, W), " & ",
    Lambda(prim_exLab_noSS, res_exLab_noSS, W), " \\\hline \end{tabular}"))
end;

=#
