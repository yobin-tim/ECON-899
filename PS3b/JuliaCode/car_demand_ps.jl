#########################################################################################
## convert from ox to julia
#########################################################################################

using DataFrames

#***************************************************************************************/
#* Main dataset + Panel structure */
#***************************************************************************************/

mPanelCharact = DataFrame(StatFiles.load("../data/Car_demand_characteristics_spec1.dta"))

n = nrow(mPanelCharact) # 6103

vYear = unique(select(mPanelCharact, [:Year])) |> Matrix # 1985-2015

T = size(vYear,1) # 31

#***************************************************************************************/
#* Outcome variables */
#***************************************************************************************/

vShare = select(mPanelCharact, [:share]) |> Matrix

vDelta_iia = select(mPanelCharact, [:delta_iia]) |> Matrix

vDelta0 = vDelta_iia

vPrice = select(car_data, [:price]) |> Matrix

