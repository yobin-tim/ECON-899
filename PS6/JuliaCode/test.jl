include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

market_clearing(prim, res)

<<<<<<< HEAD
0.95 * 7.0 + (1 - 0.95) * 15
=======







using Plots
theme(:juno)

plot(res.W_val)
scatter!(res.W_val)
plot(res.x_opt)
scatter!(res.x_opt)
>>>>>>> ps5_simulaiton_alternative
