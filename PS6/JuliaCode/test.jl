include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

market_clearing(prim, res)

res.p