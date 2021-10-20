include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

market_clearing(prim, res)

<<<<<<< HEAD
<<<<<<< HEAD
0.95 * 7.0 + (1 - 0.95) * 15
=======
=======

>>>>>>> ce279ab46af26041889ae89c207f304e1a7474cd







using Plots
theme(:juno)

plot(res.W_val)
scatter!(res.W_val)
plot(res.x_opt)
scatter!(res.x_opt)
<<<<<<< HEAD
>>>>>>> ps5_simulaiton_alternative
=======

>>>>>>> ce279ab46af26041889ae89c207f304e1a7474cd
