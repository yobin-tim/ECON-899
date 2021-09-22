<<<<<<< HEAD

#keyword-enabled structure to hold model primitives
@everywhere @with_kw struct Primitives
    β::Float64 = 0.9932 #discount rate
    α::Float64 = 1.5 #capital share
    S_grid::Array{Float64,1} = [1, 0.5] #Earnings when employed and unemployed
    Π::Array{Float64,2} = [.97 .5; .03 .5] #Transition Matrix between employment and unemployment
    na::Int64 = 700 #Number of asset grid points
    a_grid::Array{Float64,1} = collect(range(-2.0,length=na,5.0))
end

#structure that holds model results
@everywhere mutable struct Results
    val_func::Array{Float64, 2} #value function
    pol_func::Array{Float64, 2} #policy function
    q::Float64
    q_Bounds::Array{Float64,1}
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.na,2) #initial value function guess
    pol_func = zeros(prim.na,2) #zeros(prim.nk,2) #initial policy function guess
    q = (prim.β+1)/2
    q_Bounds=[prim.β, 1]
    res = Results(val_func, pol_func, q, q_Bounds) #initialize results struct
    prim, res #return deliverables
end

#Making a function for the inner loop
@everywhere module IL

    function Find_ap(S_index,a_index,res,prim)
        #unpack model primitives
            a_grid, β, α, na, Π, S_grid = prim.a_grid, prim.β, prim.α,
                prim.na, prim.Π, prim.S_grid #unpack model primitives
        #Utility Function
            function U(x)
                if x<0
                    return -Inf
                else
                    return (x^(1-α)-1)/(1-α)
                end
            end
        #Exploiting Monotonicity of V
        budget=S_grid[S_index] + a_grid[a_index];
        #Search for val_max in the found interval.
        max_val,max_ap=-Inf,0
        for ap_index=1:na
            val=U(budget - res.q*a_grid[ap_index]) +
                β*transpose(Π[:,S_index])*[res.val_func[ap_index,1]; res.val_func[ap_index,2]]
            if val>max_val
                max_val=val;
                max_ap=ap_index;
            elseif val<max_val
                break #The Value function is now declining
            end
        end
        return [max_val, max_ap]
    end
end


#Bellman Operator
function Bellman(prim::Primitives,res::Results)
    @unpack na, a_grid = prim
    v_next = zeros(na,2) #next guess of value function to fill
    for S_index=1:2
        out=pmap(a_index -> IL.Find_ap(S_index,a_index,res,prim),1:na)
        for a_index=1:na #Unpacking the pmap results
            v_next[a_index,S_index]=out[a_index][1]
            res.pol_func[a_index,S_index]=out[a_index][2]
        end
    end
    v_next #return next guess of value function
end



#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter
    while err>tol #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func))/abs(v_next[prim.na, 1]) #reset error level
        res.val_func = .8*v_next+.2*res.val_func #update value function
        n+=1
        if mod(n,50)==0
            println("Value Function iteration $(n), Error $(err)")
        end
    end
    println("Value function converged in ", n, " iterations.")
end

#Market clearing for assets/sets a new q
function MC_assets(prim,res; dist_tol::Float64 = 1e-6, dist_err::Float64 = 100.0, ES_tol=1e-2, Done=false)
    @unpack Π, na, a_grid= prim
    TransMat=zeros(2*na,2*na) #The first na points are for employed folks, and the next na are for unemployed
    for a_index=1:na
        TransMat[Int64(res.pol_func[a_index,1]),a_index]=Π[1,1] #Savings choice for those moving from emp->emp
        TransMat[Int64(res.pol_func[a_index,1])+na,a_index]=Π[2,1] #Savings choice for those moveing from emp->unemp
        TransMat[Int64(res.pol_func[a_index,2]),a_index+na]=Π[1,2] #Savings choice for those moving from unemp->emp
        TransMat[Int64(res.pol_func[a_index,2])+na,a_index+na]=Π[2,2] #Savings choice for those moving from unemp->emp
    end
    Dist=ones(2*na)*(1/(2*na)) #
    Dist_new=copy(Dist)
    while dist_err>dist_tol
        for i=1:20  #Iterate until we reach the steady-state distribution
            Dist_new=TransMat*Dist_new
        end
        dist_err=abs.(maximum(Dist_new.-Dist))
        Dist=copy(Dist_new)
    end
    #Find Excess Supply and reset q
        ExcessSupply=transpose(Dist)*vcat(a_grid,a_grid)

        if abs(ExcessSupply)>ES_tol
            #Do variant of Bisection Method
            if ExcessSupply<0
                res.q_Bounds[2]=res.q
                #Weight slightly toward old q to avoid wild fluctuations
                res.q=res.q_Bounds[1]*.3+ res.q_Bounds[2]*.7
            else
                res.q_Bounds[1]=res.q
                res.q=res.q_Bounds[1]*.7 +res.q_Bounds[2]*.3
            end
            print("Excess Supply: $(ExcessSupply), q:$(res.q)
")
        else
            Done=true
        end
        return Done
end

#solve the model
function Solve_model() #prim::Primitives, res::Results)
    prim, res = Initialize()
    converged=false
    Outer_loop_Iter=1
    while ~converged && Outer_loop_Iter<1000
        println("Beginning Asset Clearing Loop $(Outer_loop_Iter)")
        V_iterate(prim, res)
        converged=MC_assets(prim,res)
        Outer_loop_Iter+=1
    end
    return prim, res
end


#Get Distribution for Plotting
function FindDist_ForPlot(prim,res; dist_tol::Float64 = 1e-6, dist_err::Float64 = 100.0,)
    @unpack Π, na, a_grid= prim
    TransMat=zeros(2*na,2*na) #The first na points are for employed folks, and the next na are for unemployed
    for a_index=1:na
        TransMat[Int64(res.pol_func[a_index,1]),a_index]=Π[1,1] #Savings choice for those moving from emp->emp
        TransMat[Int64(res.pol_func[a_index,1])+na,a_index]=Π[2,1] #Savings choice for those moveing from emp->unemp
        TransMat[Int64(res.pol_func[a_index,2]),a_index+na]=Π[1,2] #Savings choice for those moving from unemp->emp
        TransMat[Int64(res.pol_func[a_index,2])+na,a_index+na]=Π[2,2] #Savings choice for those moving from unemp->emp
    end
    Dist=ones(2*na)*(1/(2*na)) #
    Dist_new=copy(Dist)
    while dist_err>dist_tol
        for i=1:20  #Iterate until we reach the steady-state distribution
            Dist_new=TransMat*Dist_new
        end
        dist_err=abs.(maximum(Dist_new.-Dist))
        Dist=copy(Dist_new)
    end
    return Dist, Dist[1:na].+Dist[(na+1):2*na]
end
##############################################################################
=======

#keyword-enabled structure to hold model primitives
@everywhere @with_kw struct Primitives
    β::Float64 = 0.9932 #discount rate
    α::Float64 = 1.5 #capital share
    S_grid::Array{Float64,1} = [1, 0.5] #Earnings when employed and unemployed
    Π::Array{Float64,2} = [.97 .5; .03 .5] #Transition Matrix between employment and unemployment
    na::Int64 = 700 #Number of asset grid points
    a_grid::Array{Float64,1} = collect(range(-2.0,length=na,5.0))
end

#structure that holds model results
@everywhere mutable struct Results
    val_func::Array{Float64, 2} #value function
    pol_func::Array{Float64, 2} #policy function
    q::Float64
    q_Bounds::Array{Float64,1}
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.na,2) #initial value function guess
    pol_func = zeros(prim.na,2) #zeros(prim.nk,2) #initial policy function guess
    q = (prim.β+1)/2
    q_Bounds=[prim.β, 1]
    res = Results(val_func, pol_func, q, q_Bounds) #initialize results struct
    prim, res #return deliverables
end

#Making a function for the inner loop
@everywhere module IL

    function Find_ap(S_index,a_index,res,prim)
        #unpack model primitives
            a_grid, β, α, na, Π, S_grid = prim.a_grid, prim.β, prim.α,
                prim.na, prim.Π, prim.S_grid #unpack model primitives
        #Utility Function
            function U(x)
                if x<0
                    return -Inf
                else
                    return (x^(1-α)-1)/(1-α)
                end
            end
        #Exploiting Monotonicity of V
        budget=S_grid[S_index] + a_grid[a_index];
        #Search for val_max in the found interval.
        max_val,max_ap=-Inf,0
        for ap_index=1:na
            val=U(budget - res.q*a_grid[ap_index]) +
                β*transpose(Π[:,S_index])*[res.val_func[ap_index,1]; res.val_func[ap_index,2]]
            if val>max_val
                max_val=val;
                max_ap=ap_index;
            elseif val<max_val
                break #The Value function is now declining
            end
        end
        return [max_val, max_ap]
    end
end


#Bellman Operator
function Bellman(prim::Primitives,res::Results)
    @unpack na, a_grid = prim
    v_next = zeros(na,2) #next guess of value function to fill
    for S_index=1:2
        out=pmap(a_index -> IL.Find_ap(S_index,a_index,res,prim),1:na)
        for a_index=1:na #Unpacking the pmap results
            v_next[a_index,S_index]=out[a_index][1]
            res.pol_func[a_index,S_index]=out[a_index][2]
        end
    end
    v_next #return next guess of value function
end



#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter
    while err>tol #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func))/abs(v_next[prim.na, 1]) #reset error level
        res.val_func = .8*v_next+.2*res.val_func #update value function
        n+=1
        if mod(n,50)==0
            println("Value Function iteration $(n), Error $(err)")
        end
    end
    println("Value function converged in ", n, " iterations.")
end

#Market clearing for assets/sets a new q
function MC_assets(prim,res; dist_tol::Float64 = 1e-6, dist_err::Float64 = 100.0, ES_tol=1e-2, Done=false)
    @unpack Π, na, a_grid= prim
    TransMat=zeros(2*na,2*na) #The first na points are for employed folks, and the next na are for unemployed
    for a_index=1:na
        TransMat[Int64(res.pol_func[a_index,1]),a_index]=Π[1,1] #Savings choice for those moving from emp->emp
        TransMat[Int64(res.pol_func[a_index,1])+na,a_index]=Π[2,1] #Savings choice for those moveing from emp->unemp
        TransMat[Int64(res.pol_func[a_index,2]),a_index+na]=Π[1,2] #Savings choice for those moving from unemp->emp
        TransMat[Int64(res.pol_func[a_index,2])+na,a_index+na]=Π[2,2] #Savings choice for those moving from unemp->emp
    end
    Dist=ones(2*na)*(1/(2*na)) #
    Dist_new=copy(Dist)
    while dist_err>dist_tol
        for i=1:20  #Iterate until we reach the steady-state distribution
            Dist_new=TransMat*Dist_new
        end
        dist_err=abs.(maximum(Dist_new.-Dist))
        Dist=copy(Dist_new)
    end
    #Find Excess Supply and reset q
        ExcessSupply=transpose(Dist)*vcat(a_grid,a_grid)

        if abs(ExcessSupply)>ES_tol
            #Do variant of Bisection Method
            if ExcessSupply<0
                res.q_Bounds[2]=res.q
                #Weight slightly toward old q to avoid wild fluctuations
                res.q=res.q_Bounds[1]*.3+ res.q_Bounds[2]*.7
            else
                res.q_Bounds[1]=res.q
                res.q=res.q_Bounds[1]*.7 +res.q_Bounds[2]*.3
            end
            print("Excess Supply: $(ExcessSupply), q:$(res.q)
")
        else
            Done=true
        end
        return Done
end

#solve the model
function Solve_model() #prim::Primitives, res::Results)
    prim, res = Initialize()
    converged=false
    Outer_loop_Iter=1
    while ~converged && Outer_loop_Iter<1000
        println("Beginning Asset Clearing Loop $(Outer_loop_Iter)")
        V_iterate(prim, res)
        converged=MC_assets(prim,res)
        Outer_loop_Iter+=1
    end
    return prim, res
end


#Get Distribution for Plotting
function FindDist_ForPlot(prim,res; dist_tol::Float64 = 1e-6, dist_err::Float64 = 100.0,)
    @unpack Π, na, a_grid= prim
    TransMat=zeros(2*na,2*na) #The first na points are for employed folks, and the next na are for unemployed
    for a_index=1:na
        TransMat[Int64(res.pol_func[a_index,1]),a_index]=Π[1,1] #Savings choice for those moving from emp->emp
        TransMat[Int64(res.pol_func[a_index,1])+na,a_index]=Π[2,1] #Savings choice for those moveing from emp->unemp
        TransMat[Int64(res.pol_func[a_index,2]),a_index+na]=Π[1,2] #Savings choice for those moving from unemp->emp
        TransMat[Int64(res.pol_func[a_index,2])+na,a_index+na]=Π[2,2] #Savings choice for those moving from unemp->emp
    end
    Dist=ones(2*na)*(1/(2*na)) #
    Dist_new=copy(Dist)
    while dist_err>dist_tol
        for i=1:20  #Iterate until we reach the steady-state distribution
            Dist_new=TransMat*Dist_new
        end
        dist_err=abs.(maximum(Dist_new.-Dist))
        Dist=copy(Dist_new)
    end
    return Dist, Dist[1:na].+Dist[(na+1):2*na]
end
##############################################################################
>>>>>>> 0621634e3580c37a0fd1a790ff6fee3604e7c60b
