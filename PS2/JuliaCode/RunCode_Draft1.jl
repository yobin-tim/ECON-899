#using Plots, LaTeXStrings #import the libraries we want
#Getting the Parellel Ready
    using Distributed #, SharedArrays
    #Re-initializing the workers
        rmprocs(workers())
        addprocs(6)
    @everywhere using Parameters
#Saving Details
    include("Compute_Draft1.jl")
#Solve the Model
     #initialize primitive and results structs
    @time out_primitives, out_results = Solve_model() #solve the model!
    @unpack val_func, pol_func = out_results
    @unpack a_grid = out_primitives

#Plot
Plots.plot(a_grid, val_func[:,1], title="Value Function", label="Employed")
    plot!(a_grid, val_func[:,2], label="Unemployed")
#Plots.savefig("JuliaCode/Value_Functions.png")
