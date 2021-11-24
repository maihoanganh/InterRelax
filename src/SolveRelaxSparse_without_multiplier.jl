function RelaxSparse_without_multiplier(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},s::Int64,k::Int64;assign="min",alg="MD",minimize=true,solver="Mosek",order=k,comp_opt_sol=false)
    
    println("**Interrupted relaxation based on Handelman's Positivstellensatz**")
    println("Relaxation order: k=",k)
    println("Sparsity order: s=",s)
    
   return RelaxSparse_without_multiplier1(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,s,k,assign=assign,alg=alg,minimize=minimize,solver=solver,comp_opt_sol=comp_opt_sol,order=order)
end


function RelaxSparse_without_multiplier1(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},s::Int64,d::Int64;assign="min",alg="MD",minimize=true,solver="Mosek",order=d,comp_opt_sol=false)
    
    
    I,p,lI=clique_decomp(n,m+l,[dg;dh],[[supp_f];supp_g;supp_h],order=order,alg=alg,minimize=minimize) 
    J,lJ,~=get_indcons(m,supp_g,I,p,lI,assign=assign)
    W,lW,~=get_indcons(l,supp_h,I,p,lI,assign=assign)
    
    
    
    m+=1
    
    lmon_g=[lmon_g;1]
    supp_g=[supp_g;[spzeros(UInt64,n,1)]]
    coe_g=[coe_g;[ones(Float64,1)]]
    
    dg=[dg;0]
    
    for t in 1:p
        J[t]=[J[t];m]
        lJ[t]+=1
    end
    
    
    
    
    
    println("  Number of cliques: p=", p)
    println("  Largest clique size: u=", maximum(lI))
    #println("  Possible largest block size: ", binomial(maximum(lI)+k,n))
    
    supp_f*=2
    
    
    for i in 1:m
        supp_g[i]*=2
    end
    
    for i in 1:l
        supp_h[i]*=2
    end
    #####
    
    v=Vector{Matrix{UInt64}}(undef,p)
    vmod=Vector{Matrix{UInt16}}(undef,p)
    
    
    if solver=="Mosek"
        model = Model(optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => false))
    elseif solver=="SDPT3"
        model=Model(SDPT3.Optimizer)
    elseif solver=="SDPNAL"
        model=Model(SDPNAL.Optimizer)
    elseif solver=="COSMO"
        model=Model(COSMO.Optimizer)
    else
        error("No SDP solver!!!")
    end
    
    
    
    sk=Vector{UInt64}(undef,p)
    sk_g=Vector{Vector{UInt64}}(undef,p)
    sk_h=Vector{Vector{UInt64}}(undef,p)
    
    r=1
    q=1
    maxsize=0
    
    block_G=Vector{Vector{Vector{Vector{Int64}}}}(undef,p)
    len_block_G=Vector{Vector{Vector{Int64}}}(undef,p)
    

    
    t_iter=1
    

    G=Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}}(undef, p)
    H=Vector{Vector{Vector{VariableRef}}}(undef, p)
    

    
    #####

    for t=1:p
        
        
        v[t]=get_basis(lI[t],d)
        
        sk[t]=binomial(d+lI[t],lI[t])

        sk_g[t]=Vector{UInt64}(undef,lJ[t])
        sk_h[t]=Vector{UInt64}(undef,lW[t])


        @fastmath @inbounds @simd for i in 1:lJ[t]
            sk_g[t][i]=binomial(d-dg[J[t][i]]+lI[t],lI[t])
        end

        @fastmath @inbounds @simd for i in 1:lW[t]
            sk_h[t][i]=binomial(d-dh[W[t][i]]+lI[t],lI[t])
        end
        
    end
    
    #########
    
    lsupp_U=sum(sk)      
    supp_U=spzeros(UInt64,n,lsupp_U)
    t_blo=0
    for t in 1:p
        supp_U[I[t],1+t_blo:sk[t]+t_blo]=2*v[t][:,1:sk[t]]
        t_blo+=sk[t]
    end
    
    supp_U=unique(supp_U,dims=2)
    supp_U=sortslices(supp_U,dims=2)               
    lsupp_U=size(supp_U,2)
    
    
    vec=spzeros(UInt64,n)
    
    cons=[AffExpr(0) for i=1:lsupp_U]
    
    #error("to check")
    
    for t=1:p   
        
        vmod[t]=mod.(v[t],2)
    
        r=1
        q=1

        block_G[t]=Vector{Vector{Vector{Int64}}}(undef,lJ[t])
        len_block_G[t]=Vector{Vector{Int64}}(undef,lJ[t])
        for i in 1:lJ[t]
            block_G[t][i]=Vector{Vector{Int64}}(undef,sk_g[t][i])
            len_block_G[t][i]=Vector{Int64}(undef,sk_g[t][i])
            for j in 1:sk_g[t][i]
                block_G[t][i][j]=[]
                len_block_G[t][i][j]=0
                r=j
                while len_block_G[t][i][j] <= s-1 && r <= sk_g[t][i]
                    if norm(vmod[t][:,j]-vmod[t][:,r],1)==0
                        append!(block_G[t][i][j],r)
                        len_block_G[t][i][j]+=1
                    end
                    r+=1
                end

                q=1
                while !issubset(block_G[t][i][j],block_G[t][i][q]) && q<=j-1
                    q+=1
                end

                if q<j
                    block_G[t][i][j]=[]
                    len_block_G[t][i][j]=0
                end
                #println(block_G[i][j])
                if maxsize<len_block_G[t][i][j]
                    maxsize=len_block_G[t][i][j]
                end
                #println(maxsize)
            end
        end
    

    
 
        

        G[t]=Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}(undef, lJ[t])
        H[t]=Vector{Vector{VariableRef}}(undef, lW[t])



        for i=1:lJ[t]
            G[t][i]=Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}(undef, sk_g[t][i])
            for j in 1:sk_g[t][i]


                if len_block_G[t][i][j]>=1
                    if len_block_G[t][i][j]==1
                        G[t][i][j]=@variable(model, lower_bound=0)
                        for z=1:lmon_g[J[t][i]]
                            vec=spzeros(UInt64,n)
                            vec[I[t]]=supp_g[J[t][i]][I[t],z]+2*v[t][:,block_G[t][i][j]]
                            @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-coe_g[J[t][i]][z]*G[t][i][j])
                        end
                    else 
                        G[t][i][j]=@variable(model,[1:len_block_G[t][i][j],1:len_block_G[t][i][j]],PSD)
                        for p in 1:len_block_G[t][i][j]
                            for q in p:len_block_G[t][i][j]
                                for z in 1:lmon_g[J[t][i]]
                                    if p==q
                                        vec=spzeros(UInt64,n)
                                        vec[I[t]]=v[t][:,block_G[t][i][j][p]]+v[t][:,block_G[t][i][j][q]]+supp_g[J[t][i]][I[t],z]
                                        @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-coe_g[J[t][i]][z]*G[t][i][j][p,q])
                                    else
                                        vec=spzeros(UInt64,n)
                                        vec[I[t]]=v[t][:,block_G[t][i][j][p]]+v[t][:,block_G[t][i][j][q]]+supp_g[J[t][i]][I[t],z]
                                        @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-2*coe_g[J[t][i]][z]*G[t][i][j][p,q])
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
   
    

        for i in 1:lW[t]
            H[t][i]=@variable(model, [1:sk_h[t][i]])
            for p in 1:sk_h[t][i]
                for z in 1:lmon_h[W[t][i]]
                      vec=spzeros(UInt64,n)
                      vec[I[t]]=2*v[t][:,p]+supp_h[W[t][i]][I[t],z]
                      @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-coe_h[W[t][i]][z]*H[t][i][p])
                end
            end
        end             
        
                    
    end   
                    
                    
    
                    
    @variable(model, lambda)
                    
    cons[bfind(supp_U,lsupp_U,spzeros(UInt64,n),n)]-=lambda

    for i in 1:lmon_f
        @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,supp_f[:,i],n)],coe_f[i])
    end 
                    
    
    
    @constraint(model, cons.==0)

    @objective(model, Max, lambda)
                    
    println(" Maximal matrix size:", maxsize)                
                    
    optimize!(model)

    opt_val = value(lambda)
    println("Termination status = ", termination_status(model))
    println("Optimal value = ",opt_val)
    
    opt_sol=[Vector{Float64}([]) for t in 1:p]
    if comp_opt_sol
        
        Gr=Vector{Matrix{Float64}}(undef,p)
        for t in 1:p
            println("---------------")
            println("Clique $(t):")

            Gr[t]=zeros(Float64,sk_g[t][lJ[t]],sk_g[t][lJ[t]])

            for j in 1:sk_g[t][lJ[t]]
                if len_block_G[t][lJ[t]][j]>1
                    Gr[t][block_G[t][lJ[t]][j],block_G[t][lJ[t]][j]]+=value.(G[t][lJ[t]][j])
                elseif len_block_G[t][lJ[t]][j]==1
                    Gr[t][block_G[t][lJ[t]][j],block_G[t][lJ[t]][j]]+=[value.(G[t][lJ[t]][j])]
                end

            end


            opt_sol[t]=extract_optimizer_clique(Gr[t],sk_g[t][lJ[t]],v[t][:,1:sk_g[t][lJ[t]]],lI[t],lJ[t],lW[t],lmon_g[J[t]],supp_g[J[t]],coe_g[J[t]],lmon_h[W[t]],supp_h[W[t]],coe_h[W[t]])
            println("---------------")
        end
        
    end

    return opt_val,opt_sol

end        
            
        