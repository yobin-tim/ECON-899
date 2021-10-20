include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

market_clearing(prim, res)









using Plots
theme(:juno)

plot(res.W_val)
scatter!(res.W_val)
plot(res.x_opt)
scatter!(res.x_opt)

