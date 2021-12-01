# Load packages
using LinearAlgebra, Parameters

# Auxiliary functions
# Choice probability function
function choice_probability(δ::Array{Float64}, μ::Array{Float64}; eval_jacobian::Bool = false)
    
    # number of individuals and choicesm
    r, R = size(μ);

    # Compute choice probabilities
    Λ = exp.( δ .+ μ  )
    Λ = Λ./(1 .+ sum(Λ, dims=1))
    σ = sum(Λ, dims=2)/R

    if eval_jacobian
        # Compute Jacobian
        Δ = 1/R * ((I(r) .* (σ * (1 .- σ)')) - ((1 .- I(r)) .* (σ * σ'))) ./ σ 
        return σ, Δ
    else
        return σ, nothing 
    end

end

# segment_data by market for demand estimation
function segment_data(model, market)

    # Get market id column
    market_id_col = model.market_id

    # Filter data by market
    data = model.inv_dem_est[model.inv_dem_est[!, market_id_col] .== market, :]

    # Get the observed market shares
    S = data.share
    # Get the observed prices
    P = data.price
    # Get the income levels
    Y = model.Y
    # Get the inital guess for the inverse demand
    δ = data.δ
    
    return S, P, Y, δ
end

# Model Structures
# Primitives
@with_kw struct Primitives
    λₚ_range :: Array{Float64} = [0, 1]
end


# Model
mutable struct Model 
    # Parameters
    parameters      ::Primitives                # Parameters of the model

    # Data
    market_id       :: Any                      # Market id column
    product_id      :: Any                      # Product id column
    X               :: Array{Float64, 2}        # Matrix of covariates
    Z               :: Array{Float64, 2}        # Matrix of instruments
    Y               :: Array{Float64, 2}        # Matrix of simulated data
    inv_dem_est     :: DataFrame                # DataFrame of for demand estimation
    
    # GMM estimation
    ρ               :: Array{Float64}           # GMM Residuals
end



# Demand inverter
function inverse_demand(model::Model, λₚ::Float64, market; method::String="Newton", max_iter = Inf)

    # Check the method
    valid_methods = ["Newton", "Contraction Mapping"]
    @assert (method ∈ valid_methods)

    # Get the data
    S, P, Y, δ = segment_data(model, market)

    # Compute the matrix μ[i,j] = λₚ * Y[i] * P[j]
    μ = λₚ * repeat(Y', length(S), 1) .* P

    # Compute the inverse demand
    
    # Initial guess for the inverse demand
    δ₀ = copy( δ )
    δ₁ = copy( δ )
    # Initialize the iteration counter
    iter = 0
    # Initialize the error
    err = 100
    eval_jacobian = false

    ε = 1e-12
    ε₁ = ( method == "Newton" ) ? 1 : -Inf

    # Iterate until convergence

    err_list = []
    method_flag = "Contraction Mapping"

    while (err > ε) && (iter < max_iter)
        # Compute the choice probability and the Jacobian
        if (method == "Newton") && (err < ε₁)
            eval_jacobian = true
            method_flag = "Newton"
        end

        σ, Δ = choice_probability(δ₀, μ, eval_jacobian=eval_jacobian)

        # Compute the inverse demand
        if (method == "Newton") && (err < ε₁)
            δ₁ = δ₀ + inv(Δ) * (log.(S) - log.(σ))
        else
            δ₁ = δ₀ + log.(S) - log.(σ)
        end
        
        # Update the error
        err = maximum( abs.(δ₁ - δ₀) )
        push!(err_list, err)
        # Update the inverse demand
        δ₀ = copy(δ₁)
        
        # Update the iteration counter
        iter = iter + 1
        if iter % 1000 == 0
            println("Iteration = $iter, Method = $method_flag , error = $err, tolerance = $ε, error > tolerance = $(err > ε)")     
        end

    end
    # println("Iteration = $iter, Method = $method_flag, error = $err, tolerance = $ε, error > tolerance = $(err > ε), θ = $θ")     
    market_id_col = model.market_id
    model.inv_dem_est[model.inv_dem_est[!, market_id_col] .== market, :δ] .= δ₁[:, 1]

    return err_list

end


## reference: ox code 
# /* This function evaluates the idiosyntractic component of utility */
# value(const aMu,const vParam,const t)
# {
#   decl i;
#   decl rowid=aProductID[t];
#   decl mMu=new matrix[rows(rowid)][Sim];
#   for(i=0;i<rows(vParam);i++) mMu+=vParam[i]*(aZ[t])[i];
#   aMu[0]=exp(mMu);
#   return 1;
# }
# demand(const mMu,const aShare,const aJac,const vDelta,const t,const vParam)
# {
#   decl i;
#   decl rowid=aProductID[t];
#   decl eV=exp(vDelta[rowid]).*mMu;
#   decl mS=eV./(1+sumc(eV));
#   decl vShat=meanr(mS);
#   if(aJac[0]) {
#     decl mD=diag(meanr(mS.*(1-mS)))-setdiagonal(mS*mS'/Sim,0);
#     aJac[0]=mD;
#     }
#   aShare[0]=vShat;
#   return 1;
# }
# inverse(const aDelta, const vP,const eps1,const eps)
# {
#   decl vShat,vDelta=vDelta0;
#   decl f=1000;
#   decl mJacobian=1;
#   decl rowid,t;
#   decl it,maxit=1000;
#   decl vIT=new matrix[T][1];
#   decl mMu;
#   decl time0=timer();
#   decl mNorm=<>;
#   decl maxT=T;
#   if(iprint) maxT=1;
#   parallel for(t=0;t<maxT;t++) /* Parallelize the inversion across markets. When possible use Nb processors = T (31) */
#     {
#       time0=timer();
#       /* Pre-compute the heterogeneity component of utility (indepedent of delta) */
#       value(&mMu,vP,t);       
#       rowid=aProductID[t];
#       vIT[t]=0;      
#       f=1000;
#       do{
# 	/* Evalute the demand without the jacobian matrix if the norm is larger than 1  */
# 	if(norm(f)>eps1) {
# 	  mJacobian=0;
# 	  demand(mMu,&vShat,&mJacobian,vDelta,t,vP);
# 	  f=log(vShare[rowid])-log(vShat); /* Zero function: Log difference in shares */	
# 	  vDelta[rowid]=vDelta[rowid]+f; /* Contraction mapping step */
# 	}
# 	/* Evalute the demand with the jacobian matrix if the norm is less than 1 */	
# 	else {
# 	  mJacobian=1;
# 	  demand(mMu,&vShat,&mJacobian,vDelta,t,vP);
# 	  f=log(vShare[rowid])-log(vShat); /* Zero function: Log difference in shares */
# 	  vDelta[rowid]=vDelta[rowid]+invert(mJacobian./vShat)*f; /* Newton step */
# 	}
# 	vIT[t]+=1;
# 	if(iprint==1 && t==0) {
# 	  mNorm~=(norm(f)|vIT[t]);
# 	  println("t = ",t," it ",vIT[t]," norm : ",norm(f));
# 	}
# 	//
#       }while(norm(f)>eps && vIT[t]<maxit);
#       if(norm(f)>eps) vDelta[rowid]=constant(.NaN,rowid);
#     }
#   if(iprint) {
#     DrawXMatrix(0,mNorm[0][],"Zero function norm",mNorm[1][],"Iteration");
#     ShowDrawWindow();
#     if(eps1>0) SaveDrawWindow("Inverse_iteration_newton.pdf");
#     else SaveDrawWindow("Inverse_iteration_contraction.pdf");
#   }    
#   //println("cpu time: ",(timer()-time0)/100);
#   aDelta[0]=vDelta;
#   return 1;
# }
