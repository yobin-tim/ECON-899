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
end # Primitives
@with_kw mutable struct  Results
    b       ::Array{Float64} = [0,0]
end # Primitives


function TrueData(prim)
    @unpack T, σ0, x0, ρ0 = prim
    ϵ = rand(Normal(0,sqrt(σ0)), T)
    xt=zeros(T+1)
    xt[1]=x0
    for t=2:(T+1)
        xt[t]=ρ0*xt[t-1]+ϵ[t-1]
    end
    xt=xt[2:end]
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
function J(g,W,b)
    if b[2]<0
        return 100
    else
        return transpose(g)*W*g
    end
end

function StepsAThroughD()
    prim=Primitives()
    e=eDrawsForModel(prim)
    for Exercise=4:5
        if Exercise==4
            FindM2=FindM2_MeanVar
        elseif Exercise==5
            FindM2=FindM2_VarCoVar
        end
        td=TrueData(prim)

        #Part a: Graph in three Dimensions
            @unpack ρgrid,σgrid,gpoints=prim
            Jgrid=zeros(gpoints,gpoints)
            for ρi=1:gpoints,σi=1:gpoints
                md=ModelData(prim,e, [ρgrid[ρi] σgrid[σi]])
                Jgrid[ρi,σi]=J(FindM2(td).-FindM2(md),I,[ρgrid[ρi] σgrid[σi]])
            end
            Solution=optimize(b->J(FindM2(td).-FindM2(ModelData(prim,e,b)),I,b),
                [.3 1.2],NelderMead())
            bHat1TH=Solution.minimizer
            plot(ρgrid,σgrid,Jgrid, st=:surface,
                title=L"\hat{b}^{1}_{TH}=[%$(round(bHat1TH[1],digits=4)) , %$(round(bHat1TH[2],digits=4)) ]")
            savefig("PS7\\Figures\\Exercise$(Exercise).png")
        #Part b: Find the
    end
end
