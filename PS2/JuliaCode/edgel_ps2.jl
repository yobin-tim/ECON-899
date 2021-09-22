using Parameters, Plots, LinearAlgebra, StatsPlots #import the libraries we want
# cd("C:/Users/edgel/Google Drive/UW-Madison/f21/econ899/Q1/problem_sets/PS2/JuliaCode")
include("edgel_model_functions.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res, θ = 0.5) #solve the model!
@unpack val_func, pol_func, μ, q̄ = res
@unpack a_grid = prim

##############Make plots

#stationary distribution 
Q = sum(res.μ, dims = 1)
Plots.plot(a_grid, transpose(Q))
Plots.savefig("./PS2/Figures/02_stationary_pdf_edgel.png")
Plots.plot(a_grid, transpose(μ))

cumQ = cumsum(Q, dims = 2)
Plots.plot(a_grid, transpose(cumQ))
Plots.savefig("./PS2/Figures/02_stationary_cdf_edgel.png")

#value function
Plots.plot(a_grid, transpose(val_func), title="Value Functions",
                label = ["S = e" "S = u"], legend=:topleft)
Plots.savefig("./PS2/Figures/02_Value_Functions_edgel.png")

#policy functions
Plots.plot(a_grid, transpose(pol_func), title="Policy Functions",
                label = ["S = e" "S = u"], legend=:topleft)
Plots.savefig("./PS2/Figures/02_Policy_Functions_edgel.png")

println("All done!")
################################

# I'm adding a section to get your figures publication ready
cd("/home/mitchv34/Work/2nd Year/ECON 899 (Computational Methods)/1st Quarter/Problem Sets/Shared Repo/Shared Repo") 
using LaTeXStrings
theme(:wong)
pgfplotsx()
upscale = 2 #2x upscaling in resolution
default(size=(800*upscale,600*upscale)) #Plot canvas size
fntsm = Plots.font("sans-serif", pointsize=round(10.0*upscale))
fntlg = Plots.font("sans-serif", pointsize=round(14.0*upscale))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm, tickfontsize=fntsm)

# Stationary pdf
fig_pdf = plot(a_grid, transpose(Q), xlabel = "a", ylabel = L"\mu(a)", title="Stationary Distribution", label="",
        tickfontsize = upscale * 8)

savefig(fig_pdf, "./PS2/Figures/02_pubready_stationary_pdf_edgel.pdf")

# Stationary cdf
fig_cdf = plot(a_grid, transpose(cumQ), xlabel = "a", ylabel = L"\mu(a)", title="Stationary Distribution", label="",
        tickfontsize = upscale * 8)

savefig(fig_cdf, "./PS2/Figures/02_pubready_stationary_cdf_edgel.pdf")

# Value Function
fig_val = plot(a_grid, transpose(val_func), title="Value Functions",
        xlabel="a", ylabel="v(a,s)", label = ["S = e" "S = u"], legend=:topleft, tickfontsize = upscale * 8)

savefig(fig_val, "./PS2/Figures/02_pubready_valfun_edgel.tex")

# Policy function
fig_polfun = plot(a_grid, transpose(pol_func), title="Policy Functions", xlabel="a", ylabel="g(a,s)",
        label = ["S = e" "S = u"], legend=:topleft, tickfontsize = upscale * 8)

savefig(fig_polfun, "./PS2/Figures/02_pubready_polfun_edgel.tex")        
