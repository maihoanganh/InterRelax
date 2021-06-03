
function RelaxDense(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,s::Int64)
    
    m+=1
    
    lmon_g=[lmon_g;1]
    supp_g=[supp_g;[zeros(UInt64,n,1)]]
    coe_g=[coe_g;[ones(Float64,1)]]
    
    dg=[dg;0]
    
    df=1+Int64(maximum([sum(supp_f[:,i]) for i in 1:lmon_f]))
    
    supp_f*=2
    
  
    
    lmon_thetak=Int64(1)
    supp_thetak=zeros(UInt64,n,1)
    coe_thetak=ones(Float64,1)
    
    supp_theta=2*[spzeros(UInt64,n,1) sparse(I,n,n)]
    
    for i in 1:k
        lmon_thetak,supp_thetak,coe_thetak=mulpoly(n,lmon_thetak,supp_thetak,coe_thetak,n+1,supp_theta,ones(Float64,n+1))
    end
  
    
    
    lmon_thetakf,supp_thetakf,coe_thetakf=mulpoly(n,lmon_thetak,supp_thetak,coe_thetak,lmon_f,supp_f,coe_f)
    
    v=get_basis(n,k+df)
        
    supp_U=2*v
    
    
    supp_U=sortslices(supp_U,dims=2)
    lsupp_U=size(supp_U,2)   
   
     
    sk=binomial(k+df+n,n)
    sk_g=Vector{UInt64}(undef,m)
    sk_h=Vector{UInt64}(undef,l)


    @fastmath @inbounds @simd for i in 1:m
        sk_g[i]=binomial(k+df-dg[i]+n,n)
        supp_g[i]*=2
    end
    
    @fastmath @inbounds @simd for i in 1:l
        sk_h[i]=binomial(k+df-dh[i]+n,n)
        supp_h[i]*=2
    end
    
    
   
    
    r=1
    q=1
    maxsize=0
    
    block_G=Vector{Vector{Vector{Int64}}}(undef,m)
    len_block_G=Vector{Vector{Int64}}(undef,m)
    for i in 1:m
        block_G[i]=Vector{Vector{Int64}}(undef,sk_g[i])
        len_block_G[i]=Vector{Int64}(undef,sk_g[i])
        for j in 1:sk_g[i]
            block_G[i][j]=[]
            len_block_G[i][j]=0
            r=j
            
            while len_block_G[i][j] <= s-1 && r <= sk_g[i]
                if findfirst(isodd,v[:,j]+v[:,r])==nothing
                    append!(block_G[i][j],r)
                    len_block_G[i][j]+=1
                end
                r+=1
            end
           
            q=1
            while !issubset(block_G[i][j],block_G[i][q]) && q<=j-1
                q+=1
            end
                
            if q<j
                block_G[i][j]=[]
                len_block_G[i][j]=0
            end
            #println(block_G[i][j])
            if maxsize<len_block_G[i][j]
                maxsize=len_block_G[i][j]
            end
        end
    end
        
        
   
    
    println("Maximal matrix size:", maxsize)
    
    
    #error()
        

    
    
    
    

    model=Model(with_optimizer(Mosek.Optimizer, QUIET=false))
    cons=[AffExpr(0) for i=1:lsupp_U]

    G=Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}(undef, m)
    H=Vector{Vector{VariableRef}}(undef, l)



    for i=1:m
        G[i]=Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}(undef, sk_g[i])
        for j in 1:sk_g[i]
            
            
            if len_block_G[i][j]>=1
                if len_block_G[i][j]==1
                    G[i][j]=@variable(model, lower_bound=0)
                    for z=1:lmon_g[i]
                        @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,supp_g[i][:,z]+2*v[:,block_G[i][j]],n)],coe_g[i][z]*G[i][j])
                    end
                else 
                    G[i][j]=@variable(model,[1:len_block_G[i][j],1:len_block_G[i][j]],PSD)
                    for p in 1:len_block_G[i][j]
                        for q in p:len_block_G[i][j]
                            for z in 1:lmon_g[i]
                                if p==q
                                    @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,v[:,block_G[i][j][p]]+v[:,block_G[i][j][q]]+supp_g[i][:,z],n)],coe_g[i][z]*G[i][j][p,q])
                                else
                                    @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,v[:,block_G[i][j][p]]+v[:,block_G[i][j][q]]+supp_g[i][:,z],n)],2*coe_g[i][z]*G[i][j][p,q])
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


    for i in 1:lmon_thetakf
        cons[bfind(supp_U,lsupp_U,supp_thetakf[:,i],n)]-=coe_thetakf[i]
    end
    
    @variable(model, lambda)

    for i in 1:lmon_thetak
        cons[bfind(supp_U,lsupp_U,supp_thetak[:,i],n)]+=coe_thetak[i]*lambda
    end
    
    @constraint(model, cons.==0)
    @objective(model, Max, lambda)
    optimize!(model)

    opt_val = value(lambda)
    println("Termination status = ", termination_status(model))
    println("Optimal value = ",opt_val)

    return opt_val

end        
            
         