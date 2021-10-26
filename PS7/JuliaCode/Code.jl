using Parameters, Optim, Distributions, LinearAlgebra, Plots, LaTeXStrings

@with_kw mutable struct  Primitives
    T       ::Int64             = 200        #Time Series Length
    H       ::Int64             = 10     #Number of Simulations
    ρ0      ::Float64           = .5
    σ0      ::Float64           = 1
    x0      ::Float64           = 0
    gpoints ::Int64             = 10
    ρgrid   ::Array{Float64}    = collect(range(0.35, length = gpoints, stop = 0.65))
    σgrid   ::Array{Float64}    = collect(range(0.8, length = gpoints, stop = 1.2))
    iT      ::Int64             = 4
end # Primitives
@with_kw mutable struct  Results
    bHat1TH       ::Array{Float64} = [0,0]
    bHat2TH       ::Array{Float64} = [0,0]
    e             ::Array{Float64}
    JTest         ::Float64        =100
    td            ::Array{Float64} = [0, 0] #The true data
end # Primitives


function TrueData(prim)
    @unpack T, σ0, x0, ρ0 = prim
    ϵ = rand(Normal(0,sqrt(σ0)), T)
    xt=zeros(T+1)
    xt[1]=x0
    for t=2:(T+1)
        xt[t]=ρ0*xt[t-1]+ϵ[t-1]
    end
    # xt=xt[2:end]
    return xt
end
function eDrawsForModel(prim)
    @unpack H,T=prim
    e=zeros(T,H)
    for hi=1:H
        e[:,hi]=rand(Normal(0,1),T)
    end
    return e
end
function ModelData(prim,e,b)
    @unpack T, σ0, x0, ρ0, H = prim
    if b[2]<0
        #Here J is set to return a very high value so this will not be picked
        return zeros(T,H)
    else
        yt=zeros(T+1,H)
        for hi=1:H
            for t=2:(T+1)
                yt[t,hi]=b[1]*yt[t-1,hi]+sqrt(b[2])*e[t-1,hi]
            end
        end
        yt=yt[2:end,:]
        return yt
    end
end
function FindM2_MeanVar(data)
    M=zeros(2)
    if size(data,2)>1 #We are dealing with the Model Data
        Hmeans=sum(data,dims=1)./size(data,1)
        M[1]=sum(Hmeans)/size(data,2)
        #data.-Hmeans does the subtraction for each row
        Hvars=sum((data.-Hmeans).^2,dims=1)./size(data,1)
        M[2]=sum(Hvars)/size(data,2)
    else
        M[1]=sum(data)/length(data)
        M[2]=sum((data.-M[1]).^2)/length(data)
    end
    return M
end
function m_MeanVar(x,ind, mean)
    return [x[ind]  (x[ind]-mean)^2]
end
function FindM2_VarCoVar(data)
    M=zeros(2)
    if size(data,2)>1 #We are dealing with the Model Data
        Hmeans=sum(data,dims=1)./size(data,1)
        #data.-Hmeans does the subtraction for each row
        Hvars=sum((data.-Hmeans).^2,dims=1)./size(data,1)
        M[1]=sum(Hvars)/size(data,2)
        HCovars=sum((data[2:end,:].-Hmeans).*(data[1:(end-1),:].-Hmeans),dims=1)./(size(data,1)-1)
        M[2]=sum(HCovars)/size(data,2)
    else #We are dealing with the true data
        Mean=sum(data)/length(data)
        M[1]=sum((data.-Mean).^2)/length(data)
        M[2]=sum((data[2:end].-Mean).*(data[1:(end-1)].-Mean))/(length(data)-1)
    end
    return M
end
function m_VarCovar(x,ind,mean)
    return [(x[ind]-mean)^2      (x[ind]-mean)*(x[ind-1]-mean)]
end
function FindM3(data)
    M=zeros(3)
    if size(data,2)>1 #We are dealing with the Model Data
        Hmeans=sum(data,dims=1)./size(data,1)
        M[1]=sum(Hmeans)/size(data,2)
        #data.-Hmeans does the subtraction for each row
        Hvars=sum((data.-Hmeans).^2,dims=1)./size(data,1)
        M[2]=sum(Hvars)/size(data,2)
        HCovars=sum((data[2:end,:].-Hmeans).*(data[1:(end-1),:].-Hmeans),dims=1)./(size(data,1)-1)
        M[3]=sum(HCovars)/size(data,2)
    else #We are dealing with the true data
        Mean=sum(data)/length(data)
        M[1]=Mean;
        M[2]=sum((data.-Mean).^2)/length(data)
        M[3]=sum((data[2:end].-Mean).*(data[1:(end-1)].-Mean))/(length(data)-1)
    end
    return M
end
function m3(x,ind,mean)
    [x[ind]   (x[ind]-mean)^2      (x[ind]-mean)*(x[ind-1]-mean)]
end
function J(g,W,b) #The objective function
    if b[2]<0 #We shouldn't allow there to be a negative variance term
        return 1e25
    else
        return transpose(g)*W*g
    end
end

function GraphAndFindbHat(W,prim,res,FindM,Exercise; NeweyWest=false)
    @unpack ρgrid,σgrid,gpoints=prim
    @unpack e, td = res
    Jgrid=zeros(gpoints,gpoints)
    for ρi=1:gpoints,σi=1:gpoints
        md=ModelData(prim,e, [ρgrid[ρi] σgrid[σi]])
        Jgrid[ρi,σi]=J(FindM(td).-FindM(md),W,[ρgrid[ρi] σgrid[σi]])
    end
    Solution=optimize(b->J(FindM(td).-FindM(ModelData(prim,e,b)),W,b),
        [.3 1.2],NelderMead())
    bHat=Solution.minimizer
    plot(ρgrid,σgrid,Jgrid, st=:surface,
        title=L"\hat{b}^{1}_{TH}=[%$(round(bHat[1],digits=4)) , %$(round(bHat[2],digits=4)) ]")
    if NeweyWest
        # savefig("PS7\\Figures\\Exercise$(Exercise)NeweyWestCorrection.png")
        savefig("PS7/Figures/Exercise$(Exercise)NeweyWestCorrection.png")
    else
        # savefig("PS7\\Figures\\Exercise$(Exercise).png")
        savefig("PS7/Figures/Exercise$(Exercise).png")
    end
    return bHat
end

function NeweyWest(prim::Primitives,simdata::Array{Float64},
        m::Function, MTH::Array{Float64})
    @unpack iT,H,T = prim
    Hmeans=sum(simdata,dims=1)./size(simdata,1)
    #Defining Γ function
        function ΓjTH(j)
            out=zeros(length(MTH),length(MTH))
            for hi=1:H, ti=(j+1+1):T
                #Necessary to add one twice above because one moment includes a lag
                out+=(m(simdata[:,hi],ti,Hmeans[hi]).-MTH)*
                    transpose(m(simdata[:,hi],ti-j,Hmeans[hi]).-MTH)
            end
            return (1/((T-1)*H))*out
        end
    #Finding STH
        SyTH=ΓjTH(0)
        for ji=1:iT
            SyTH+=(1-ji/(iT+1))*(ΓjTH(ji)+transpose(ΓjTH(ji)))
        end
        STH=(1+1/H)*SyTH
    #Returning WStar
        return inv(STH)
end

function Find∇g(res,prim, FindM; s=1e-7)
    @unpack e, bHat2TH = res
    ∂g∂ρ=(FindM(ModelData(prim,e,bHat2TH)).-FindM(ModelData(prim,e,bHat2TH-[s 0])))./s
    ∂g∂σ=(FindM(ModelData(prim,e,bHat2TH)).-FindM(ModelData(prim,e,bHat2TH-[0 s])))./s
    return hcat(∂g∂ρ,∂g∂σ)
end

function StepsAThroughD()
    prim=Primitives()
    res=Results(e=eDrawsForModel(prim))
    for Exercise=4:6
        if Exercise==4
            FindM=FindM2_MeanVar
            m=m_MeanVar
        elseif Exercise==5
            FindM=FindM2_VarCoVar
            m=m_VarCovar
        elseif Exercise==6
            FindM=FindM3
            m=m3
        end
        res.td=TrueData(prim)
        #Function for parts a and b

        #Part a: Graph in three Dimensions
            res.bHat1TH=GraphAndFindbHat(I,prim,res,FindM,Exercise)
            print("
________________________________________________________________________________
                    Results for Exercise $(Exercise)
________________________________________________________________________________\n")
            println("The estimate of b using W = I is ", res.bHat1TH,".")

        #Part b: Use NeweyWest to update your guess of bHat
            md_bHat1TH=ModelData(prim,res.e, res.bHat1TH)
            WStar=NeweyWest(prim,ModelData(prim,res.e, md_bHat1TH),
                    m, FindM(md_bHat1TH))
            res.bHat2TH=GraphAndFindbHat(WStar,prim,res,FindM,Exercise,
                NeweyWest=true)
            println("The estimate of b using Wstar is ", res.bHat2TH,".")

        #Part c
            ∇g=Find∇g(res,prim,FindM)
                print("∇g= \n\n")
                display(∇g)
            VarCovarbHat2TH=(1/prim.T)*inv(transpose(∇g)*WStar*∇g)
                print("\n The variance-covariance matrix for bHat2TH is given by \n\n")
                display(VarCovarbHat2TH)
            StdErrorsbHat2TH=sqrt.(diag(VarCovarbHat2TH))
                print("\n The Standard errors are given by \n\n")
                display(StdErrorsbHat2TH)
        #Part d: Computing the Value of the J test
            res.JTest=prim.T*(prim.H/(1+prim.H))*
                J(FindM(res.td).-FindM(ModelData(prim,res.e,res.bHat2TH)),
                    WStar,res.bHat2TH)
                print("\n The J-Test is $(res.JTest)")
    end #Exercise Loop
end #End Function StepsAThroughD
