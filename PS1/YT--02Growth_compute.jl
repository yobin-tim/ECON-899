@time begin
using Parameters, Plots #import the libraries we want
include("YT--02Growth_model.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func = res
@unpack k_grid = prim

##############Make plots
#value function
Plots.plot(k_grid, val_func, title="Value Function", labels = ["Z_g" "Z_b"])
Plots.savefig("02_Value_Functions.png")

#policy functions
Plots.plot(k_grid, pol_func, title="Policy Functions", labels = ["Z_g" "Z_b"])
Plots.savefig("02_Policy_Functions.png")

#changes in policy function
pol_func_δ = copy(pol_func).-k_grid
Plots.plot(k_grid, pol_func_δ, title="Policy Functions Changes", labels = ["Z_g" "Z_b"])
Plots.savefig("02_Policy_Functions_Changes.png")

println("All done!")
################################
end
