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
    @unpack a_grid, na, S_grid = out_primitives

#Plotting results
using Plots, LaTeXStrings #import the libraries we want
Plots.plot(a_grid, val_func[:,1], title="Value Function", label="Employed")
    plot!(a_grid, val_func[:,2], label="Unemployed")
    Plots.savefig("Value_Functions.png")
    #Plotting Policy functions
        function PolicyPolots()
            a_hat=0
            for ai=1:na
                if a_grid[Int64.(pol_func[ai,1])]<=a_grid[ai]
                    a_hat=[a_grid[ai]];
                    break
                end
            end
            Plots.plot(a_grid, a_grid[Int64.(pol_func[:,1])], title="Policy Functions", label="Employed")
                plot!(a_grid, a_grid[Int64.(pol_func[:,2])], label="Unemployed")
                plot!(a_grid,a_grid, label="45 line", legend=:bottomright)
                vline!(a_hat, label=L"\hat{a}",color=:purple)
                annotate!(a_hat[1]-.15, 4.5, text("$(round(a_hat[1],digits=3))", :purple, :right, 12))
                Plots.savefig("Policy_Functions.png")
        end
        PolicyPolots()
    #Plotting Distribution
        function DistPlots()
            TS_Distribution, SS_WealthDistribution=FindDist_ForPlot(out_primitives,out_results)
            MaxNonZero=1
            ForDistPlot=copy(TS_Distribution)
            for i=1:na
                if ForDistPlot[i]==0
                    ForDistPlot[i]=NaN
                end
                if ForDistPlot[na+i]==0
                    ForDistPlot[na+i]=NaN
                end
                if SS_WealthDistribution[i]!=0
                    MaxNonZero=copy(i)
                end
            end
            Plots.plot(a_grid[1:MaxNonZero], ForDistPlot[1:MaxNonZero], title="Distribution of Assets
where q=$(round(out_results.q,digits=8))",
                label="Employed")
                plot!(a_grid[1:MaxNonZero], ForDistPlot[(na+1):na+MaxNonZero], label="Unemployed", xlabel="Assets")
                Plots.savefig("Distribution.png")
            #Lorenz Curve
            n_lorenz=1000
            Lorenz=zeros(n_lorenz,2)
                Lorenz[:,1]=collect(range(0,length=n_lorenz,1)) #First column is percent of population
                i=1
                for a_index=1:na
                    if sum(SS_WealthDistribution[1:a_index])<=Lorenz[i,1]
                        Lorenz[i,2]=Lorenz[i,2]+TS_Distribution[a_index]*(a_grid[a_index]+S_grid[1]) +
                            TS_Distribution[na+a_index]*(a_grid[a_index]+S_grid[2]) #Second column is cumulative assets
                    else
                        while sum(SS_WealthDistribution[1:a_index])>Lorenz[i,1]
                            i+=1
                            Lorenz[i,2]=Lorenz[i-1,2]+0; #copy over the previous cumulative wealth
                        end
                    end
                    Lorenz[i,2]=Lorenz[i,2]+TS_Distribution[a_index]*(a_grid[a_index]+S_grid[1]) +
                        TS_Distribution[na+a_index]*(a_grid[a_index]+S_grid[2])
                end
                #Calculating Gini
                Gini=sum(Lorenz[:,1].-Lorenz[:,2])/(sum(Lorenz[:,1].-Lorenz[:,2])+sum(Lorenz[:,1]))
                #Lorenz[:,2]=Lorenz[:,2]./Lorenz[n_lorenz,2] #express cumulative assets as a percentage
                #print(Lorenz)
                Plots.plot(100*Lorenz[:,1],100*Lorenz[:,2], title="Lorenz Curve.
The Gini Coefficient is $(round(Gini,digits=8))",
                xlabel="% of Population",
                    ylabel="% of Assets", legend=:bottomright, label="Lorenz")
                    plot!(100*Lorenz[:,1],100*Lorenz[:,1], label="Line of Equality")
                    Plots.savefig("Lorenz.png")
        end
        DistPlots()
#
