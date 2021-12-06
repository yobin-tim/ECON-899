mutable struct Results

    est :: Array{Float64,1}

end

function load_data()

    mPanelCharact = DataFrame(StatFiles.load("../data/Car_demand_characteristics_spec1.dta"))

    n = nrow(mPanelCharact) # 6103

    vYear = sort(unique(mPanelCharact[!, :Year]))

    T = size(vYear,1) # 31

    vShare = select(mPanelCharact, [:share]) |> Matrix

    vDelta = select(mPanelCharact, [:delta_iia]) |> Matrix

    vPrice = select(mPanelCharact, [:price]) |> Matrix
 
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
 
    #* Load price and differentation IVs *
  
    mDemandIV = DataFrame(StatFiles.load("../data/Car_demand_iv_spec1.dta"))

    aIVlist = [:i_import,
               :diffiv_local_0,
               :diffiv_local_1,
               :diffiv_local_2,
               :diffiv_local_3,
               :diffiv_ed_0]

    mExclIV = mDemandIV[!,aIVlist] |> Matrix

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

    mEta = DataFrame(StatFiles.load("../data/Simulated_type_distribution.dta")) |> Matrix

    return vYear, vShare, vDelta, aProductID, mX, mZ, mEta, mIV
    
end

function Initialize()

    est = zeros(6103)

    res = Results(est)

    return res
    
end

function value(vParam, t)

    aZ = mZ[aProductID[t]] .* mEta' 

    aMu = zeros(size(aProductID[t],1), 100)

    aMu = exp.(vParam * aZ)

    return aMu

end

function demand(aMu, aJac, vDelta, t)

    rowid = aProductID[t]

    eV = exp.(vDelta[rowid]) .* aMu

    tmp = (1 .+ sum(eV, dims = 1))

    mS = eV./ tmp

    vShat = mean(mS, dims = 2)

    if aJac == 1
        
        mD = diagm(0 => vec(mean(mS .* (1 .- mS), dims = 2))) -
            mS*mS'/100 -diagm(0 => diag(mS*mS'/100))

    else

        mD = 0

    end        

    return vShat, mD
    
end

function inverse(res::Results, aDelta, vParam)

    @unpack est = res

    for t = 1:31
        
        mMu = value(vParam, t)

        rowid = aProductID[t]

        err_list = []

        err = 100

        iter = 0
        
        tmp = 0

        converge = 0

        f = 100

        while converge == 0 

            if err > 1

                tmp = demand(mMu, 0, aDelta, t)

                vShat = tmp[1];

                f = log.(vShare[rowid]) - log.(vShat);

                aDelta[rowid] = aDelta[rowid] + f;

            else
                
                tmp = demand(mMu, 1, aDelta, t)

                vShat = tmp[1];
                
                mJacobian = tmp[2];

                f = log.(vShare[rowid]) - log.(vShat);

                aDelta[rowid] = aDelta[rowid] + inv(mJacobian./vShat)*f;
 
            end

            err = norm(f)
            
            iter += 1
            
            if err < 10^(-12)

                res.est[rowid] = aDelta[rowid]

                converge = 1

            end

        end

    end

    return res
    
end

function ivreg(mY, mVariables, mInstruments, mWeight)

    mInvXZX = inv((mVariables'*mInstruments)*mWeight*(mInstruments'*mVariables))

    vIV=mInvXZX*(mVariables'mInstruments)*mWeight*(mInstruments'mY);  

    return vIV
    
end

function gmm_obj(vParam)

    res = Initialize();
    
    vDelta_init = vDelta;
    
    inverse(res, vDelta_init, vParam)

    vLParam = ivreg(res.est, mX, mIV, A)

    vXi = res.est - mX * vLParam

    mG = vXi .* mIV

    g = sum(mG, dims = 1)

    adFunc = g*A*g'/100

    return adFunc[1]
    
end

