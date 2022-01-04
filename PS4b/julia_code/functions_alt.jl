using Parameters, DelimitedFiles, LinearAlgebra,  DataFrames, Optim

@with_kw mutable struct Primitives
    # Parameters of the model
    α     ::Int64       = 2
    λ     ::Int64       = -4      #Stockout penalty
    β     ::Float64     = .99
    euler ::Float64     = 2.71828   #Euler's number
    # State space and Transition Matrix
    S   ::Array{Float64,2}  = convert(Array{Float64,2},
        readdlm("ox_code/PS4_state_space.csv", ',', skipstart=1)[:,3:end]) #I, c, p
        NoStates ::Int64 = size(S,1)
    Fa0 ::Array{Float64,2} =convert(Array{Float64,2},
        readdlm("ox_code/PS4_transition_a0.csv", ',', skipstart=1)[:,3:end])
    Fa1 ::Array{Float64,2} =convert(Array{Float64,2},
        readdlm("ox_code/PS4_transition_a1.csv", ',', skipstart=1)[:,3:end])
    #Functions
    F   ::Function = (a::Int64) -> (a==1) ? Fa1 : Fa0
end


@with_kw mutable struct Results
    EV_P1 :: Array{Float64,1}
    EV_CCP :: Array{Float64,1}
    EV_Findingλ :: Array{Float64,1}
    λOpt :: Float64
end

function InitializeModel()
    prim=Primitives()
    res=Results(EV_P1=zeros(prim.NoStates),
        EV_CCP=zeros(prim.NoStates),
        EV_Findingλ=zeros(prim.NoStates),
        λOpt=prim.λ)
    return prim, res
end

function Bellman(prim::Primitives,res::Results; tol=1e-4, FindingOptimalλ::Bool=false,
        Altλ::Float64=-4.0)
    @unpack F, euler, β, NoStates, S,α,λ =prim
    @unpack EV_P1=res
    if FindingOptimalλ #For problem 4
        λ=copy(Altλ)
    end
    function U_NoShock(a::Int64,si::Int64) #si is the state index
        if a==1
            return α*S[si,2]-S[si,3]
        elseif a==0 && S[si,1]>0
            return α*S[si,2]
        elseif a==0 && S[si,1]==0
            return λ*S[si,2]
        else
            print("Error")
        end
    end
    function EV_Iter(EV::Array{Float64,1})
        EV_inner=zeros(NoStates)
        for si=1:NoStates
            EV_inner[si]=log(sum(exp.(
                    [U_NoShock(ai,si)+β*dot(F(ai)[si,:],EV) for ai=0:1]
                    )))+euler
        end
        return EV_inner
    end
    err,iter,EV_old,EV=100,1,EV_P1,zeros(NoStates)
    while err>tol
        EV=EV_Iter(EV_old)
        err=maximum(abs.(EV-EV_old))
        EV_old=copy(EV)
        if iter%100==0 && ~FindingOptimalλ
            println("Value function Iteration $(iter)")
        end
        iter+=1
    end
    if FindingOptimalλ
        return EV, U_NoShock
    else
        return EV
    end
end

function PrepareData()
    Data=convert(DataFrame,
        readdlm("ox_code/PS4_simdata.csv", ',', skipstart=1)[:,2:3]) #Choice, State Index
        rename!(Data, [:choice,:state_id])
        Data[!,:choice] = convert.(Float64,Data[!,:choice])
        Data[!,:state_id] = convert.(Int64,Data[!,:state_id])
    return Data
end


function Estimate_CCP(prim::Primitives,Data::DataFrame)
    @unpack NoStates =prim

    CCP=zeros(NoStates)
    for si=1:NoStates
        Sample=filter(row->row.state_id==si,Data)
        CCP[si]=sum(Sample[:,:choice]) / maximum( [nrow(Sample) 1.0] )
        #Constrain to lie between 0.001 and 0.999
            if CCP[si]<0.001
                CCP[si]=0.001
            elseif CCP[si]>.999
                CCP[si]=0.999
            end
    end
    return CCP
end

function FindEV_CCP(prim::Primitives,Data::DataFrame)
    @unpack F,β,NoStates,euler,S,λ,α=prim
    function U_NoShock(a::Int64,si::Int64) #si is the state index
        if a==1
            return α*S[si,2]-S[si,3]
        elseif a==0 && S[si,1]>0
            return α*S[si,2]
        elseif a==0 && S[si,1]==0
            return λ*S[si,2]
        else
            print("Error")
        end
    end
    function U_vec(a::Int64)
        return [U_NoShock(a,si) for si=1:NoStates]
    end
    CCP=Estimate_CCP(prim,Data)
    e_vec ::Function =(a::Int64) -> (a==1) ? euler.-log.(CCP) : euler.-log.((1 .-CCP))
    FP=F(0).*(1 .- CCP)+F(1).*CCP
    #Id=Matrix{Int}(I, NoStates, NoStates)
    EV_CCP=inv(I - β*FP) *   ((1 .- CCP).*(U_vec(0)+e_vec(0))+
                                CCP.*(U_vec(1)+e_vec(1)))
    return EV_CCP
end

function LogLiklihoodOfλ(λ::Float64,prim::Primitives,res::Results, Data::DataFrame)
    @unpack NoStates,β,F = prim
    res.EV_Findingλ, U_NoShock = Bellman(prim,res, FindingOptimalλ=true, Altλ=λ)
    function Prob_a1()
        v_tilde=zeros(NoStates)
        for si=1:NoStates
            v_tilde[si]=(U_NoShock(1,si)+β*dot(F(1)[si,:],res.EV_Findingλ))-
                        (U_NoShock(0,si)+β*dot(F(0)[si,:],res.EV_Findingλ))
        end
        return (1 .+ exp.(-v_tilde)).^(-1)
    end
    liklihood=0.0
    P=Prob_a1()
    for si=1:NoStates
        Sample=filter(row->row.state_id==si,Data)
        if nrow(Sample)>0
            for ri=1:nrow(Sample)
                liklihood+=Sample[ri,:choice]*log(P[si])+(1-Sample[ri,:choice])*log(1-P[si])
            end
        end
    end
    return liklihood
end

function Findingλ(prim::Primitives,res::Results, Data::DataFrame)
    λOpt=optimize(λInner -> -LogLiklihoodOfλ(λInner[1] , prim,res, Data),
        [Float64(prim.λ)],method=BFGS()).minimizer
    return λOpt[1]
end

function SolveModel()
    prim,res=InitializeModel()
    #Problem 1
        res.EV_P1=Bellman(prim,res)
    #Problem 2
        Data=PrepareData()
        res.EV_CCP=FindEV_CCP(prim,Data)
    #Problem 4
        res.λOpt=Findingλ(prim,res,Data)
    return res
end
