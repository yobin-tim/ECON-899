#########################################################################################
## convert from ox to julia
#########################################################################################

using StatFiles, DataFrames, Statistics, LinearAlgebra

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

## TODO: save different length of vector in array##
aProductID = []
for i = 1:T
    aProductID[i] = mPanelCharact[(mPanelCharact.Year .== vYear[i]), :].index
end
###################################################

aProductID = mPanelCharact[(mPanelCharact.Year .== 1985.0), :].index




#/***************************************************************************************/
# /* Random coefficients */
#/***************************************************************************************/

mEta = DataFrame(StatFiles.load("../data/Simulated_type_distribution.dta")) |> Matrix

Sim = nrow(mEta)

## TODO:different size of matrix

aZ = mZ[aProductID] .* mEta' # 149 × 100 for 1985



