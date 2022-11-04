
function RelaxDense_without_multiplier(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,s::Int64;solver="Mosek",L=1.0)
    
    println("**Interrupted relaxation based on Handelman's Positivstellensatz**")
    
    println("Relaxation order: k=",k)
    println("Factor width upper bound: s=",s)
    
    m+=1
    
    lmon_g=[lmon_g;1]
    supp_g=[supp_g;[zeros(UInt64,n,1)]]
    coe_g=[coe_g;[ones(Float64,1)]]
    
    dg=[dg;0]
    
    
    supp_f*=2
    
  
    lmon_bcons_power=Vector{Int64}(undef,k+1)
    supp_bcons_power=Vector{Matrix{UInt64}}(undef,k+1)
    coe_bcons_power=Vector{Vector{Float64}}(undef,k+1)
    
    lmon_bcons_power[1]=Int64(1)
    supp_bcons_power[1]=zeros(UInt64,n,1)
    coe_bcons_power[1]=ones(Float64,1)
    
    lmon_bcons=n+1
    supp_bcons=[spzeros(UInt64,n,1) 2*sparse(I,n,n)]
    coe_bcons=[L;-ones(Float64,n)]
    
    
    for i in 1:k
        lmon_bcons_power[i+1],supp_bcons_power[i+1],coe_bcons_power[i+1]=mulpoly(n,lmon_bcons_power[i],supp_bcons_power[i],coe_bcons_power[i],lmon_bcons,supp_bcons,coe_bcons)
    end
    
    
    
    v=get_basis(n,k)
        
    supp_U=2*v
    
    
    supp_U=sortslices(supp_U,dims=2)
    lsupp_U=size(supp_U,2)   
   
     
    sk=binomial(k+n,n)
    sk_g=Vector{Vector{Int64}}(undef,m)
    sk_h=Vector{Int64}(undef,l)

    for i in 1:m
        sk_g[i]=Vector{Int64}(undef,k-dg[i]+1)
        for r in 0:k-dg[i]
            sk_g[i][r+1]=binomial(k-dg[i]-r+n,n)
        end
        supp_g[i]*=2
    end
    
    
    @fastmath @inbounds @simd for i in 1:l
        sk_h[i]=binomial(k-dh[i]+n,n)
        supp_h[i]*=2
    end
    
    
   vmod=mod.(v,2)
    
    r=1
    q=1
    maxsize=0
    
    block_G=Vector{Vector{Vector{Vector{Int64}}}}(undef,m)
    len_block_G=Vector{Vector{Vector{Int64}}}(undef,m)
    for i in 1:m
        block_G[i]=Vector{Vector{Vector{Int64}}}(undef,k-dg[i]+1)
        len_block_G[i]=Vector{Vector{Int64}}(undef,k-dg[i]+1)
        for t in 0:k-dg[i]
            block_G[i][t+1]=Vector{Vector{Int64}}(undef,sk_g[i][t+1])
            len_block_G[i][t+1]=Vector{Int64}(undef,sk_g[i][t+1])
            for j in 1:sk_g[i][t+1]
                block_G[i][t+1][j]=[]
                len_block_G[i][t+1][j]=0
                r=j

                while len_block_G[i][t+1][j] <= s-1 && r <= sk_g[i][t+1]
                    #if all(el->iseven(el)==true, v[:,j]+v[:,r])#
                    if norm(vmod[:,j]-vmod[:,r],1)==0
                        append!(block_G[i][t+1][j],r)
                        len_block_G[i][t+1][j]+=1
                    end
                    r+=1
                end

                q=1
                while !issubset(block_G[i][t+1][j],block_G[i][t+1][q]) && q<=j-1
                    q+=1
                end

                if q<j
                    block_G[i][t+1][j]=[]
                    len_block_G[i][t+1][j]=0
                end
                #println(block_G[i][j])
                if maxsize<len_block_G[i][t+1][j]
                    maxsize=len_block_G[i][t+1][j]
                end
            end
        end
    end
        
        
   
    
    println("Maximal matrix size:", maxsize)
    
    
    #error()
        
    #ENV["MATLAB_ROOT"] = "/usr/local/MATLAB/R2018a/toolbox/local"
    
    if solver=="Mosek"
        model=Model(optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => false))
    elseif solver=="SDPT3"
        model=Model(SDPT3.Optimizer)
    elseif solver=="SDPNAL"
        model=Model(SDPNAL.Optimizer)
    elseif solver=="COSMO"
        model=Model(COSMO.Optimizer)
    else
        error("No SDP solver!!!")
    end
    
    
    cons=[AffExpr(0) for i=1:lsupp_U]

    G=Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}}(undef, m)
    H=Vector{Vector{VariableRef}}(undef, l)



    for i=1:m
        G[i]=Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}(undef, k-dg[i]+1)
        for r in 0:k-dg[i]
            G[i][r+1]=Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}(undef, sk_g[i][r+1])
            for j in 1:sk_g[i][r+1]
                if len_block_G[i][r+1][j]>=1
                    if len_block_G[i][r+1][j]==1
                        G[i][r+1][j]=@variable(model, lower_bound=0)
                        for z=1:lmon_g[i]
                            for a=1:lmon_bcons_power[r+1]
                                @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,supp_g[i][:,z]+supp_bcons_power[r+1][:,a]+2*v[:,block_G[i][r+1][j]],n)],coe_g[i][z]*coe_bcons_power[r+1][a]*G[i][r+1][j])
                            end
                        end
                    else 
                        G[i][r+1][j]=@variable(model,[1:len_block_G[i][r+1][j],1:len_block_G[i][r+1][j]],PSD)
                        for p in 1:len_block_G[i][r+1][j]
                            for q in p:len_block_G[i][r+1][j]
                                for z in 1:lmon_g[i]
                                    for a=1:lmon_bcons_power[r+1]
                                        if p==q
                                            @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,v[:,block_G[i][r+1][j][p]]+v[:,block_G[i][r+1][j][q]]+supp_g[i][:,z]+supp_bcons_power[r+1][:,a],n)],coe_g[i][z]*G[i][r+1][j][p,q]*coe_bcons_power[r+1][a])
                                        else
                                            @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,v[:,block_G[i][r+1][j][p]]+v[:,block_G[i][r+1][j][q]]+supp_g[i][:,z]+supp_bcons_power[r+1][:,a],n)],2*coe_g[i][z]*G[i][r+1][j][p,q]*coe_bcons_power[r+1][a])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
   
    

    for i in 1:l
        H[i]=@variable(model, [1:sk_h[i]])
        for p in 1:sk_h[i]
            for z in 1:lmon_h[i]
                  @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,2*v[:,p]+supp_h[i][:,z],n)],coe_h[i][z]*H[i][p])
            end
        end
    end


    for i in 1:lmon_f
        cons[bfind(supp_U,lsupp_U,supp_f[:,i],n)]-=coe_f[i]
    end
    
    @variable(model, lambda)

    cons[bfind(supp_U,lsupp_U,zeros(UInt64,n),n)]+=lambda
    
    @constraint(model, cons.==0)
    @objective(model, Max, lambda)
    optimize!(model)

    opt_val = value(lambda)
    println("Termination status = ", termination_status(model))
    println("Primal status = ", primal_status(model))
    println("Optimal value = ",opt_val)

    return opt_val

end        
            
         