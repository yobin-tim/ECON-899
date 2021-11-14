# Loading Packages
using Parameters, LinearAlgebra, Plots, Latexify, DataFrames, LaTeXStrings

# Loadifng Programs 
include("./hopenhayn_rogerson.jl")

# Set plot theme
theme(:vibrant) 
default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style

# Initialize the model's parameters and results struct
prim, res = Initialize();

# Create the structure for experiments 1 and 2 with random disturbances
_, res_1 = Initialize();
_, res_2 = Initialize();
_, res_3 = Initialize();

# First we solve for the case with no random disturbances
solve_model_no_dist(prim, res)

# and then we add random disturbances to the model
# Fist we will create a dictionary of that will index the result structure with the random disturbance
results = Dict(0.0 => res, 1.0 => res_1, 2.0 => res_2, 3.0 => res_3 ) # α = 0.0 means no disturbances
# then we iterate over the random disturbances and solve the model
for (α, res_struct) in results
    if α  != 0.0
        find_equilibrium(prim, res_struct, α)
    end
end


# Plot the results
p2 = plot(prim.s_vals, res.x_opt, size=(800,600),
            title="Decision Rules of Exit", label="Benchmark Model",
            linestyle =:dash, markershape = :auto, legend = :right)
xticks!(prim.s_vals)
yticks!(0:0.1:1)
xlabel!("Firm Productivity")
ylabel!("Exit Probability")


for (α, res_struct) in results
    if α  != 0.0
        plot!(prim.s_vals, res_struct.x_opt, size=(800,600),
                title="Decision Rules of Exit", label="TV1 Shocks α = $α ",
                linestyle =:dash, markershape = :auto)
    end
end
current()

# savefig(p2, "./PS6/Document/Figures/decision_rules_2.pdf")
savefig(p2, "../Document/Figures/decision_rules_2.pdf")

# Save results to a table
## Error for i 
n_opt = [prim.n_optim.(prim.s_vals, r.p) for (_, r) in results]
n_incumbents = [ prim.n_optim.(prim.s_vals, r.p)' * r.μ for (_, r) in results]
n_entrants = [ prim.n_optim.(prim.s_vals, r.p)' * prim.ν * results[i].M for (_, r) in results]
n_total = [n_incumbents[i] + n_entrants[i] for i in 1:length(n_incumbents) ]
fraction_labor_entrants = [n_entrants[i] ./ n_total[i] for i in 1:length(n_incumbents) ]

df = DataFrame(["Price Level" => [r.p for (_, r) in results]
                "Mass of Incumbents" => [sum(r.μ) for (_, r) in results]
                "Mass of Entrants" => [sum(r.M .* prim.ν) for (_, r) in results]
                "Mass of Exits" => [sum(r.μ .* r.x_opt) for (_, r) in results]
                "Aggregate Labor" => n_total
                "Labor of Incumbents" => n_incumbents
                "Labor of Entrants" => n_entrants
                "Fraction of Labor Entrants" => fraction_labor_entrants]
                )
                df = round.(df, digits=3)
                colnames = names(df)
                df[!, :id] = [( α == 0.0 ) ?  "Standard" : "TV1 Shock α = $α" for (α, _) in results]
df = unstack(stack(df, colnames), :variable, :id, :value)

laetx_table = latexify(df; env=:table, latex=false)

open("./PS6/Document/Tables/table_1.tex", "w") do file
    write(file, laetx_table)
#open("../Document/Tables/table_2.tex", "w") do file
#    write(file, laetx_table)

end
