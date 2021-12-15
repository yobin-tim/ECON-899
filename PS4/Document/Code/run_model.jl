## Changeing current directory
## cd(expanduser("~/pathtoEcon-899/ECON-899/"))
cd(expanduser("~/Box/Econ899/Git/ECON-899/")) # e.g. Hiroaki

using Distributed, SharedArrays, JLD

#add processes
workers()
addprocs(2)

@Distributed.everywhere include("./PS4/JuliaCode/conesa_kueger.jl");
#@Distributed.everywhere include("./conesa_kueger.jl");

#Exercise 1: (Problem Set 4)
@time out_K_path,out_Ft,out_vf_trans= TransitionPath(TrySaveMethod=false,Experiment=1)
#Exercise 2: (Problem Set 4)
@time out_K_path_Exp2,out_Ft_Exp2,out_vf_trans_Exp2=
    TransitionPath(TrySaveMethod=false,Experiment=2)

using Plots, LaTeXStrings

theme(:juno)
plot(1:81,out_K_path[:], ylabel="Aggregate Capital",
    xlabel="Time",label="Excercise 1")
    plot!(1:81,out_K_path_Exp2[:], legend=:bottomright, ylabel="Aggregate Capital",
        xlabel="Time",label="Excercise 2")
    savefig("./PS4/Document/Figures/ComparingTransitions.png")

function Exercise1Prob2(kpath; α=.36,δ=.06,N_final=66,J_R=46, TransitionNo=81)
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


#Exercise 3 Problem 1
out_primStart, out_resStart= MarketClearing(use_Fortran=false, tol = 1e-3);
#out_primEnd, out_resEnd= MarketClearing(ss=false, use_Fortran=false, tol = 1e-3);
function SolveProblem3()
    @unpack nZ, nA, N_final, σ = out_primStart
    γ=1
    EV=zeros(nZ,nA,N_final)
    EVj=zeros(N_final)
    PortionWhoAreMadeBetterOff=0;
    for zi=1:nZ
        for ai=1:nA
            for ji=1:N_final
                EV[zi,ai,ji]=(out_vf_trans[ai,zi,ji,1] / out_resStart.val_fun[ai,zi,ji])^(1/(γ*(1-σ)))
                EVj[ji]+=out_resEnd.F[ai,zi,ji]*EV[zi,ai,ji]
                if EV[zi,ai,ji]>1
                    PortionWhoAreMadeBetterOff+=out_resStart.F[ai,zi,ji]
                end
            end #age loop
        end #asset loop
    end #z loop
    print(EVj)
    plot(1:N_final,EVj, ylabel="Consumption Equivalent Variation",
        xlabel="Age",
        title="$(round(PortionWhoAreMadeBetterOff,digits=2)*100)% of the Population would
Support the Reform", legend=false)
        savefig("./PS4/Document/Figures/Exercise1Problem3.png")
end
SolveProblem3()
