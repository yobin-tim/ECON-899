include("./hopenhayn_rogerson.jl")

prim, res = Initialize()

market_clearing(prim, res)

0.95 * 7.0 + (1 - 0.95) * 15