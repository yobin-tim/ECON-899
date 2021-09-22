<<<<<<< HEAD
#########################################
# Problem Set 2, ECON 899
# Goal Due Date: 9/20/2021
#########################################

#ReadMe: This file calls functions from PS2 Model. Refer to that file to modify functions.

@time begin
using Parameters, Plots, LinearAlgebra #import the libraries we want
cd("/Volumes/GoogleDrive/My Drive/Courses/Taken Courses/Fall 2021/ECON 899/PS2/JuliaCode/")
include("YT-PS2 Model.jl") #import the functions that solve our growth model
prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func, μ = res
@unpack A_grid = prim

## Make plots
#value function
Plots.plot(A_grid, transpose(val_func), title="Value Functions", labels = ["s=e" "s=u"])
Plots.savefig("YT-PS2_Value_Functions.png")

#wealth distribution
Plots.plot(A_grid, transpose(pol_func), title = "Policy Functions", labels = ["s=e" "s=u"])
Plots.plot!(A_grid, A_grid, labels = "45 deg", ls = :dash)
Plots.savefig("YT-PS2_Policy_Functions.png")
end
=======
#########################################
# Problem Set 2, ECON 899
# Goal Due Date: 9/20/2021
#########################################

#ReadMe: This file calls functions from PS2 Model. Refer to that file to modify functions.

@time begin
using Parameters, Plots, LinearAlgebra #import the libraries we want
cd("/Volumes/GoogleDrive/My Drive/Courses/Taken Courses/Fall 2021/ECON 899/PS2/JuliaCode/")
include("YT-PS2 Model.jl") #import the functions that solve our growth model
prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func, μ = res
@unpack A_grid = prim

## Make plots
#value function
Plots.plot(A_grid, transpose(val_func), title="Value Functions", labels = ["s=e" "s=u"])
Plots.savefig("YT-PS2_Value_Functions.png")

#wealth distribution
Plots.plot(A_grid, transpose(pol_func), title = "Policy Functions", labels = ["s=e" "s=u"])
Plots.plot!(A_grid, A_grid, labels = "45 deg", ls = :dash)
Plots.savefig("YT-PS2_Policy_Functions.png")
end
>>>>>>> 0621634e3580c37a0fd1a790ff6fee3604e7c60b
