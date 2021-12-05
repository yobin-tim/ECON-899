#########################################################################################
## convert from ox to julia
#########################################################################################

using StatFiles, DataFrames, Statistics, LinearAlgebra, Plots

#***************************************************************************************/
#* Main dataset + Panel structure */
#***************************************************************************************/
mPanelCharact = DataFrame(StatFiles.load("../data/Car_demand_characteristics_spec1.dta"))

n = nrow(mPanelCharact) # 6103

vYear = sort(unique(mPanelCharact[!, :Year]))

T = size(vYear,1) # 31

#***************************************************************************************/
#* Outcome variables */
#***************************************************************************************/

vShare = select(mPanelCharact, [:share]) |> Matrix

vDelta_iia = select(mPanelCharact, [:delta_iia]) |> Matrix

vDelta0 = vDelta_iia

vPrice = select(mPanelCharact, [:price]) |> Matrix

#***************************************************************************************/
#* Characteristics and instrumentsal variables */
#***************************************************************************************/
#* Load characteristics */
 
varlist = [:price,
           :dpm,
           :hp2wt,
           :size,
           :turbo,
           :trans,
           :Year_1986,
           :Year_1987,
           :Year_1988,
           :Year_1989,
           :Year_1990,
           :Year_1991,
           :Year_1992,
           :Year_1993,
           :Year_1994,
	   :Year_1995,
           :Year_1996,
           :Year_1997,
           :Year_1998,
           :Year_1999,
           :Year_2000,
           :Year_2001,
           :Year_2002,
           :Year_2003,
           :Year_2004,
           :Year_2005,
	   :Year_2006,
           :Year_2007,
           :Year_2008,
           :Year_2009,
           :Year_2010,
           :Year_2011,
           :Year_2012,
           :Year_2013,
           :Year_2014,
           :Year_2015,
           :model_class_2,
	   :model_class_3,
           :model_class_4,
           :model_class_5,
           :cyl_2,
           :cyl_4,
           :cyl_6,
           :cyl_8,
           :drive_2,
           :drive_3,
           :Intercept]

exo_varlist=[:dpm,
             :hp2wt,
             :size,
             :turbo,
             :trans,
             :Year_1986,
             :Year_1987,
             :Year_1988,
             :Year_1989,
             :Year_1990,
             :Year_1991,
             :Year_1992,
             :Year_1993,
             :Year_1994,
             :Year_1995,
             :Year_1996,
             :Year_1997,
             :Year_1998,
             :Year_1999,
             :Year_2000,
             :Year_2001,
             :Year_2002,
             :Year_2003,
             :Year_2004,
             :Year_2005,
             :Year_2006,
             :Year_2007,
             :Year_2008,
             :Year_2009,
             :Year_2010,
             :Year_2011,
             :Year_2012,
             :Year_2013,
             :Year_2014,
             :Year_2015,
             :model_class_2,
             :model_class_3,
             :model_class_4,
             :model_class_5,
             :cyl_2,
             :cyl_4,
             :cyl_6,
             :cyl_8,
             :drive_2,
             :drive_3,
             :Intercept]

mX = mPanelCharact[!, varlist] |> Matrix
 
mean.(eachcol(mX))

#* Load price and differentation IVs *
  
mDemandIV = DataFrame(StatFiles.load("../data/Car_demand_iv_spec1.dta"))

aIVlist = [:i_import,
           :diffiv_local_0,
           :diffiv_local_1,
           :diffiv_local_2,
           :diffiv_local_3,
           :diffiv_ed_0]

mExclIV = mDemandIV[!,aIVlist] |> Matrix

mean.(eachcol(mExclIV))

## "~": horizontal concatention 
mIV1 = mPanelCharact[!, exo_varlist] |> Matrix

mIV = hcat(mIV1, mExclIV)

#* Non-linear attributes */
mZ = mPanelCharact[!, [:price]] |>Matrix

#/* Pre-compute the row IDs for each market */
index = rownumber.(eachrow(mPanelCharact))

insertcols!(mPanelCharact, :index => index)

#ref: https://stackoverflow.com/questions/56366632/julia-array-with-different-sized-vectors
aProductID = Array{Any}(undef,31);

for i = 1:T
    aProductID[i] = mPanelCharact[(mPanelCharact.Year .== vYear[i]), :].index
end

#/***************************************************************************************/
# /* Random coefficients */
#/***************************************************************************************/

mEta = DataFrame(StatFiles.load("../data/Simulated_type_distribution.dta")) |> Matrix

vParam = 0.6

function value(vParam, t)

    aZ = mZ[aProductID[t]] .* mEta' 

    aMu = zeros(size(aProductID[t],1), 100)

    aMu = exp.(vParam * aZ)

    return aMu

end

function demand(mMu, aJac, vDelta, t)

    rowid = aProductID[t]

    eV = exp.(vDelta[rowid]) .* mMu

    tmp = (1 .+ sum(eV, dims = 1))

    mS = eV./ tmp

    vShat = mean(mS, dims = 2)

    if aJac == 1
        
        mD = diagm(0 => vec(mean(mS .* (1 .- mS), dims = 2))) -
            mS*mS'/Sim -diagm(0 => diag(mS*mS'/Sim))

    else

        mD = 0

    end        

    return vShat, mD
    
end

## Q1-1: Contraction mapping

t = 1

vDelta = vDelta_iia

aJac = 0

eps0 = 1e-12

err_list = []

err = 100

iter = 0

while err > eps0
    
    tmp = demand(mMu, aJac, vDelta, t)

    vShat = tmp[1];

    f = log.(vShare[rowid]) - log.(vShat);

    vDelta[rowid]= vDelta[rowid]+f;

    err = norm(f)
    
    push!(err_list, err)

    iter += 1

end

x = 1:iter;

plot(x, err_list,
     title = "Contraction Mapping",
     xlabel = "Iteration",
     ylabel = "Norm",
     color =:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/Q1_contraction.pdf")


## Q1-2 Combination of contraction mapping and Newton
vDelta = vDelta_iia

function inverse(aDelta, vParam)

    for t = 1:T
        
        mMu = value(vParam, t)

        rowid = aProductID[t]

        err_list = []

        err = 100

        iter = 0
        
        tmp = 0

        while err > 10^(-12)

            if err > 1
                
                tmp = demand(mMu, 0, vDelta, t)

                vShat = tmp[1];

                f = log.(vShare[rowid]) - log.(vShat);

                vDelta[rowid] = vDelta[rowid]+f;

            else
                
                tmp = demand(mMu, 1, vDelta, t)

                vShat = tmp[1];
                
                mJacobian = tmp[2];

                f = log.(vShare[rowid]) - log.(vShat);

                vDelta[rowid] = vDelta[rowid]+inv(mJacobian./vShat)*f;

            end

            err = norm(f)
            
            push!(err_list, err)

            iter += 1
            
            ## TODO: save vDelta properly
         
        end

        println(t)

    end

end

    
x = 1:iter;

plot(x, err_list,
     title = "Contraction Mapping and Newton",
     xlabel = "Iteration",
     ylabel = "Norm",
     color =:black,
     legend = false,
     lw = 2)

savefig("../Document/Figures/Q1_combination.pdf")

## Q2-Grid serch
