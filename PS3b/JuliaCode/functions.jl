# Load packages




# Demand inversion Contraction Mapping (Barry et al., 1995)
function inverse_demand_cm(δ_0, s)
    
end

# Demand inverter
function inverse_demand(params::Array{Float64} ; method::String="Newton")

    if method=="Newton"

    elseif method=="CM"

    else
        error("Method not implemented, implemented methods are \"Newton\" and \"CM\" ")
    end

end


inverse_demand([1.0]; method="CM")