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




