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
    @unpack a_grid, na = out_primitives

#Plot
using Plots, LaTeXStrings #import the libraries we want
Plots.plot(a_grid, val_func[:,1], title="Value Function", label="Employed")
    plot!(a_grid, val_func[:,2], label="Unemployed")
    Plots.savefig("Value_Functions.png")
Plots.plot(a_grid, pol_func[:,1], title="Policy Functions", label="Employed")
    plot!(a_grid, pol_func[:,2], label="Unemployed")
    Plots.savefig("Policy_Functions.png")
Type_Specific_Distribution, SS_Distribution=FindDist_ForPlot(out_primitives,out_results)
    Plots.plot(a_grid, Type_Specific_Distribution[1:na], title="Value Function", label="Employed")
        plot!(a_grid, Type_Specific_Distribution[(na+1):2*na], label="Unemployed")
        Plots.savefig("Distribution.png")
#
