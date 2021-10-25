using Distributed, SharedArrays, JLD

#add processes
workers()
addprocs(1)

@Distributed.everywhere include("./conesa_kueger_tmp.jl");

prim_check, res_check = Convergence()

## Exercise 1 2.
using Plots, LaTeXStrings

time = [1:1:res_check.N;]

plot(time, res_check.r_path,
     title = "Interest Rate",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/interest_rate.pdf")

plot(time, res_check.w_path,
     title = "Wage Path",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/wage.pdf")

plot(time, res_check.K_path,
     title = "Capital Path",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/capital.pdf")

plot(time, res_check.L_path,
     title = "Labor Path",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/labor.pdf")

# Exercise 1.3
d = load("../Data/Initial_Conditions.jld") # Initial Conditions from PS3

val_SS_0 = d["V_0"] # Value function from the initial steady state 

val_0 = res_check.val_fun_path[:,:,:,1] #Value function from generations 

EV = (val_0./val_SS_0).^(1/(prim_check.γ*(1-prim_check.σ)))

A = EV.*d["Γ_0"]

EV_j = dropdims(sum(dropdims(sum(A, dims =2), dims =2), dims=1), dims = 1)

age_length = [1:1:66;]

plot(age_length, EV_j,
     title = "EV",
     xlabel = "Age",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/EV.pdf")

# Exercise 2
prim_anticipate, res_anticipate = Convergence(false)

time_length = [1:1:res_anticipate.N;]

plot(time_length, res_anticipate.r_path,
     title = "Interest Rate",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/interest_rate_ex2.pdf")

plot(time_length, res_anticipate.w_path,
     title = "Wage Path",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/wage_ex2.pdf")

plot(time_length, res_anticipate.K_path,
     title = "Capital Path",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/capital_ex2.pdf")

plot(time_length, res_anticipate.L_path,
     title = "Labor Path",
     xlabel = "Time",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/labor_ex2.pdf")

# Exercise 1.3

val_SS_0 = d["V_0"] # Value function from the initial steady state 

val_0 = res_anticipate.val_fun_path[:,:,:,1] #Value function from generations 

EV = (val_0./val_SS_0).^(1/(prim_check.γ*(1-prim_check.σ)))

A = EV.*d["Γ_0"]

EV_j = dropdims(sum(dropdims(sum(A, dims =2), dims =2), dims=1), dims = 1)

age_length = [1:1:66;]

plot(age_length, EV_j,
     title = "EV",
     xlabel = "Age",
     ylabel = "",
     color=:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/EV_ex2.pdf")
