function get_blocks_Cliq_mix(k::Int64,n::Int64,m::Int64,l::Int64,Usupp::SparseMatrixCSC{UInt64,Int64},lmon_g::Vector{UInt64},lmon_h::Vector{UInt64},supp_g#=::Vector{SparseMatrixCSC{UInt64}}=#,supp_h#=::Vector{SparseMatrixCSC{UInt64}}=#,coe_g::Vector{Vector{Float64}},coe_h::Vector{Vector{Float64}},v::Matrix{UInt64},sk::UInt64,sk_g::Vector{UInt64},sk_h::Vector{UInt64})
    
    #println(typeof(Usupp))
    
    Usupp=sortslices(Usupp,dims=2)
    Usupp=unique(Usupp,dims=2)
    lUsupp=size(Usupp,2)

    block_g=Vector{Vector{Vector{UInt64}}}(undef,m)
    block_h=Vector{Vector{Vector{UInt64}}}(undef,l)

    lblock_g=Vector{UInt64}(undef,m)
    lblock_h=Vector{UInt64}(undef,l)

    lt_block_g=Vector{Vector{UInt64}}(undef,m)
    lt_block_h=Vector{Vector{UInt64}}(undef,l)
    
    graph=SimpleGraph(sk)
    for p in 1:sk, q in 1:p
        if bfind(Usupp,lUsupp,v[:,p]+v[:,q],n)!=0
           add_edge!(graph,p,q)
        end
    end

    y=1
    for i in 1:m
        graph=SimpleGraph(sk_g[i])
        for p in 1:sk_g[i]
            for q in 1:p
                while y<=lmon_g[i]
                    if bfind(Usupp,lUsupp,v[:,p]+v[:,q]+supp_g[i][:,y],n)!=0
                        break
                    else
                        y+=1
                    end
                end
                if y<=lmon_g[i]
                   add_edge!(graph,p,q)
                end
                y=1
            end
        end
        block_g[i]=connected_components(graph)
        lblock_g[i]=length(block_g[i])
        lt_block_g[i]=[length(block_g[i][j]) for j in 1:lblock_g[i]]
    end



    for i in 1:l
        graph=SimpleGraph(sk_h[i])
        for p in 1:sk_h[i]
            for q in 1:p
                while y<=lmon_h[i]
                    if bfind(Usupp,lUsupp,v[:,p]+v[:,q]+supp_h[i][:,y],n)!=0
                        break
                    else
                        y+=1
                    end
                end
                if y<=lmon_h[i]
                   add_edge!(graph,p,q)
                end
                y=1
            end
        end
        block_h[i]=connected_components(graph)
        lblock_h[i]=length(block_h[i])
        lt_block_h[i]=[length(block_h[i][j]) for j in 1:lblock_h[i]]
    end
    
    #=
   for i in 1:m
        for j in 1:lblock_g[i]
            for p in 1:lt_block_g[i][j], q in 1:p
                for z in 1:lmon_g[i]
                    Usupp=[Usupp v[:,block_g[i][j][p]]+v[:,block_g[i][j][q]]+supp_g[i][:,z]]
                end
            end
        end
    end
    for i in 1:l
        for j in 1:lblock_h[i]
            for p in 1:lt_block_h[i][j], q in 1:p
                for z in 1:lmon_h[i]
                    Usupp=[Usupp v[:,block_h[i][j][p]]+v[:,block_h[i][j][q]]+supp_h[i][:,z]]
                end
            end
        end
    end
      =#
    
    Usupp=unique(Usupp,dims=2)
    lUsupp=size(Usupp,2)

    
    return Usupp,lUsupp,block_g,block_h,lblock_g,lblock_h,lt_block_g,lt_block_h

end


function run_get_blocks_Cliq_mix(n::Int64,p::Int64,t::Int64,k::Int64,lI::Vector{Int64},lJ::Vector{Int64},lW::Vector{Int64},I::Vector{Vector{UInt16}},J::Vector{Vector{UInt64}},W::Vector{Vector{UInt64}},Supp::SparseMatrixCSC{UInt64,Int64},IndA::Vector{Vector{UInt64}},lmon_g::Vector{UInt64},lmon_h::Vector{UInt64},supp_g::Array{SparseMatrixCSC{UInt64,Ti} where Ti<:Integer,1},supp_h::Array{SparseMatrixCSC{UInt64,Ti} where Ti<:Integer,1},coe_g::Vector{Vector{Float64}},coe_h::Vector{Vector{Float64}},dg::Vector{Int64},dh::Vector{Int64})
    
    
    Usupp=Vector{Matrix{UInt64}}(undef,p)
    lUsupp=Vector{UInt64}(undef,p)
    
    block_g=Vector{Vector{Vector{Vector{UInt64}}}}(undef,p)
    block_h=Vector{Vector{Vector{Vector{UInt64}}}}(undef,p)

    lblock_g=Vector{Vector{UInt64}}(undef,p)
    lblock_h=Vector{Vector{UInt64}}(undef,p)

    lt_block_g=Vector{Vector{Vector{Int64}}}(undef,p)
    lt_block_h=Vector{Vector{Vector{Int64}}}(undef,p)
    
    v=Vector{Matrix{UInt64}}(undef,p)
    sk=Vector{UInt64}(undef,p)
    sk_g=Vector{Vector{UInt64}}(undef,p)
    sk_h=Vector{Vector{UInt64}}(undef,p)
    
    
    
    for j in 1:p
        #v[j],sk[j],sk_g[j],sk_h[j]=initial_mix(lI[j],lJ[j],lW[j],k,dg[J[j]],dh[W[j]])
         

        v[j]=get_basis(lI[j],k)
        
        sk[j]=binomial(k+lI[j],lI[j])

        sk_g[j]=Vector{UInt64}(undef,lJ[j])
        sk_h[j]=Vector{UInt64}(undef,lW[j])


        @fastmath @inbounds @simd for i in 1:lJ[j]
            sk_g[j][i]=binomial(k-dg[J[j][i]]+lI[j],lI[j])
        end

        @fastmath @inbounds @simd for i in 1:lW[j]
            sk_h[j][i]=binomial(k-dh[W[j][i]]+lI[j],lI[j])
        end 
        Usupp[j],lUsupp[j],block_g[j],block_h[j],lblock_g[j],lblock_h[j],lt_block_g[j],lt_block_h[j]=get_blocks_Cliq_mix(k,lI[j],lJ[j],lW[j],[Supp[I[j],IndA[j]] 2*v[j]],lmon_g[J[j]],lmon_h[W[j]],[supp_g[i][I[j],:] for i in J[j]],[supp_h[i][I[j],:] for i in W[j]],coe_g[J[j]],coe_h[W[j]],v[j],sk[j],sk_g[j],sk_h[j])
        
        
         #Usupp[j],lUsupp[j],block_g[j],block_h[j],lblock_g[j],lblock_h[j],lt_block_g[j],lt_block_h[j]=get_blocks_Cliq_mix(k,lI[j],lJ[j],lW[j],Supp[I[j],IndA[j]],lmon_g[J[j]],lmon_h[W[j]],[supp_g[i][I[j],:] for i in J[j]],[supp_h[i][I[j],:] for i in W[j]],coe_g[J[j]],coe_h[W[j]],v[j],sk[j],sk_g[j],sk_h[j])
        
        Usupp[j]=sortslices(Usupp[j],dims=2)
        Usupp[j]=unique(Usupp[j],dims=2)
        lUsupp[j]=size(Usupp[j],2)
    end
    
    
    lSupp=sum(lUsupp)
    Supp=spzeros(UInt64,n,lSupp)
    t_blo=0
    for j in 1:p
        Supp[I[j],1+t_blo:lUsupp[j]+t_blo]=Usupp[j]
        t_blo+=lUsupp[j]
    end
    Supp=unique(Supp,dims=2)
    lSupp=size(Supp,2)
    
    IndA=[UInt64[] for i in 1:p]
    ind=zeros(UInt64,1)
    for j in 1:lSupp
        ind=findall(y->y>0,Supp[:,j])
        for i in 1:p
            if issubset(ind,I[i])
                push!(IndA[i],j)
            end
        end
    end
    
    
    
    for iter in 2:t
    
        for j in 1:p

            #println("Clique ",j,"th: ==================")
            Usupp[j],lUsupp[j],block_g[j],block_h[j],lblock_g[j],lblock_h[j],lt_block_g[j],lt_block_h[j]=get_blocks_Cliq_mix(k,lI[j],lJ[j],lW[j],Supp[I[j],IndA[j]],lmon_g[J[j]],lmon_h[W[j]],[supp_g[i][I[j],:] for i in J[j]],[supp_h[i][I[j],:] for i in W[j]],coe_g[J[j]],coe_h[W[j]],v[j],sk[j],sk_g[j],sk_h[j])
            if iter==t
                Usupp[j]=sortslices(Usupp[j],dims=2)
                Usupp[j]=unique(Usupp[j],dims=2)
                lUsupp[j]=size(Usupp[j],2)
            end
        end
        
        if iter<t
            
            lSupp=sum(lUsupp)
            Supp=spzeros(UInt64,n,lSupp)
            t_blo=0
            for j in 1:p
                Supp[I[j],1+t_blo:lUsupp[j]+t_blo]=Usupp[j]
                t_blo+=lUsupp[j]
            end
        
            Supp=unique(Supp,dims=2)
            lSupp=size(Supp,2)

            IndA=[UInt64[] for i in 1:p]
            for j in 1:lSupp
                ind=findall(y->y>0,Supp[:,j])
                for i in 1:p
                    if issubset(ind,I[i])
                        push!(IndA[i],j)
                    end
                end
            end
    #println("==================")
      
        end
    end
    
    return v,Supp,lSupp,block_g,block_h,lblock_g,lblock_h,lt_block_g,lt_block_h,sk,sk_g,sk_h
end





function RelaxSparseCSTS(n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{SparseMatrixCSC{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{SparseMatrixCSC{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::SparseMatrixCSC{UInt64},coe_f::Vector{Float64},dg::Vector{Int64},dh::Vector{Int64},k::Int64,t::Int64,s::Int64;assign="min",alg="MD",minimize=true,solver="Mosek",order=k,comp_opt_sol=false,)
    
    
    println("**Interrupted relaxation based on Handelman's Positivstellensatz**")
    println("Relaxation order: k=",k)
    println("Term sparsity order: t=",t)
    println("Upper bound on maximal matrix size: s=",s)
    
    
    
    I,p,lI=clique_decomp(n,m+l,[dg;dh],[[supp_f];supp_g;supp_h],order=order,alg=alg,minimize=minimize) 
    J,lJ,~=get_indcons(m,supp_g,I,p,lI,assign=assign)
    W,lW,~=get_indcons(l,supp_h,I,p,lI,assign=assign)
    
    m+=1
    
    lmon_g=[lmon_g;1]
    supp_g=[supp_g;[spzeros(UInt64,n,1)]]
    coe_g=[coe_g;[ones(Float64,1)]]
    
    dg=[dg;0]
    
    for y in 1:p
        J[y]=[J[y];m]
        lJ[y]+=1
    end
    
    supp_f*=2
    
    
    for i in 1:m
        supp_g[i]*=2
    end
    
    for i in 1:l
        supp_h[i]*=2
    end
    
    #dg=2*dg
    #dh=2*dh
    
    println("  Number of cliques: p=", p)
    println("  Largest clique size: u=", maximum(lI))
    
    
    
    
    lSupp=lmon_f+sum(lmon_g)+sum(lmon_h)
    Supp=spzeros(UInt64,n,lSupp)
    Supp[:,1:lmon_f]=supp_f
    t_blo=lmon_f
    for j in 1:m
        Supp[:,1+t_blo:lmon_g[j]+t_blo]=supp_g[j]
        t_blo+=lmon_g[j]
    end
    for j in 1:l
        Supp[:,1+t_blo:lmon_h[j]+t_blo]=supp_h[j]
        t_blo+=lmon_h[j]
    end
    
    #Need to remove
    IndA=[UInt64[] for i in 1:p]
    for j in 1:lSupp
        ind=findall(y->y>0,Supp[:,j])
        for i in 1:p
            if issubset(ind,I[i])
                push!(IndA[i],j)
            end
        end
    end
    
    
    v,supp_U,lsupp_U,block_g,block_h,lblock_g,lblock_h,lt_block_g,lt_block_h,sk,sk_g,sk_h=run_get_blocks_Cliq_mix(n,p,t,k,lI,lJ,lW,I,J,W,Supp,IndA,lmon_g,lmon_h,supp_g,supp_h,coe_g,coe_h,dg,dh)
        
    println("  Even symmetry: ",all(iseven.(supp_U)))
   
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
    
    
    
    
    #mex,lt_mex,l_mex,mex_cliq,lmex_cliq=get_mex_mix(n,lsupp_U,supp_U,p,I,lI)
    
    #Indf,lIndf=decomp_obj(supp_f,lmon_f,I,p,n)
    
    
    G=Vector{Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}}}(undef, p)
    H=Vector{Vector{Vector{VariableRef}}}(undef, p)
    
    supp_U_h=Vector{Vector{Matrix{UInt64}}}(undef,p)
    lsupp_U_h=Vector{Vector{Int64}}(undef,p)
    
    vec=spzeros(UInt64,n)
    
    cons=[AffExpr(0) for i=1:lsupp_U]
    
    t_iter=1
    r=1
    q=1
    
    block_G=Vector{Vector{Vector{Vector{Vector{Int64}}}}}(undef,p)
    len_block_G=Vector{Vector{Vector{Vector{Int64}}}}(undef,p)
    
    maxsize=0
    
    for y in 1:p

    
        r=1
        q=1

        block_G[y]=Vector{Vector{Vector{Vector{Int64}}}}(undef,lJ[y])
        len_block_G[y]=Vector{Vector{Vector{Int64}}}(undef,lJ[y])
        for i in 1:lJ[y]
            block_G[y][i]=Vector{Vector{Vector{Int64}}}(undef,lblock_g[y][i])
            len_block_G[y][i]=Vector{Vector{Int64}}(undef,lblock_g[y][i])
            for j in 1:lblock_g[y][i]
                block_G[y][i][j]=Vector{Vector{Int64}}(undef,lt_block_g[y][i][j])
                len_block_G[y][i][j]=Vector{Int64}(undef,lt_block_g[y][i][j])
                for a in 1:lt_block_g[y][i][j]
                    block_G[y][i][j][a]=[]
                    len_block_G[y][i][j][a]=0
                    r=a
                    while len_block_G[y][i][j][a] <= s-1 && r <= lt_block_g[y][i][j]
                        append!(block_G[y][i][j][a],block_g[y][i][j][r])
                        len_block_G[y][i][j][a]+=1
                        r+=1
                    end

                    q=1
                    while !issubset(block_G[y][i][j][a],block_G[y][i][j][q]) && q<=a-1
                        q+=1
                    end

                    if q<a
                        block_G[y][i][j][a]=[]
                        len_block_G[y][i][j][a]=0
                    end
                    #println(block_G[y][i][j][a])
                     if maxsize<len_block_G[y][i][j][a]
                         maxsize=len_block_G[y][i][j][a]
                     end
                    #println(maxsize)
                end
            end
        end
        
        
        
        
        G[y]=Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}}(undef, lJ[y])
        for i in 1:lJ[y]
            G[y][i]=Vector{Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}}(undef, lblock_g[y][i])
            for j in 1:lblock_g[y][i]
                G[y][i][j]=Vector{Union{VariableRef,Symmetric{VariableRef,Array{VariableRef,2}}}}(undef, lt_block_g[y][i][j])
                for b in 1:lt_block_g[y][i][j]
                    if len_block_G[y][i][j][b]>0
                        if len_block_G[y][i][j][b]==1
                            G[y][i][j][b]=@variable(model, lower_bound=0)
                            for z=1:lmon_g[J[y][i]]
                                vec=spzeros(UInt64,n)
                                vec[I[y]]=supp_g[J[y][i]][I[y],z]+2*v[y][:,block_G[y][i][j][b]]
                                @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-coe_g[J[y][i]][z]*G[y][i][j][b])
                            end
                        else
                            G[y][i][j][b]=@variable(model,[1:len_block_G[y][i][j][b],1:len_block_G[y][i][j][b]],PSD)
                            for a in 1:len_block_G[y][i][j][b]
                                for q in a:len_block_G[y][i][j][b]
                                    for z in 1:lmon_g[J[y][i]]
                                        if a==q
                                            vec=spzeros(UInt64,n)
                                            vec[I[y]]=v[y][:,block_G[y][i][j][b][a]]+v[y][:,block_G[y][i][j][b][q]]+supp_g[J[y][i]][I[y],z]
                                            @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-coe_g[J[y][i]][z]*G[y][i][j][b][a,q])
                                        else
                                            vec=spzeros(UInt64,n)
                                            vec[I[y]]=v[y][:,block_G[y][i][j][b][a]]+v[y][:,block_G[y][i][j][b][q]]+supp_g[J[y][i]][I[y],z]
                                            @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-2*coe_g[J[y][i]][z]*G[y][i][j][b][a,q])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        
       
        supp_U_h[y]=Vector{Matrix{UInt64}}(undef,lW[y])
        lsupp_U_h[y]=zeros(Int64,lW[y])
        
        H[y]=Vector{Vector{VariableRef}}(undef, lW[y])
        
        for i in 1:lW[y]
            supp_U_h[y][i]=zeros(UInt64,lI[y],Int(0.5*sum(lt_block_h[y][i][j]*(lt_block_h[y][i][j]+1) for j in 1:lblock_h[y][i])))
            t_iter=1
            for j in 1:lblock_h[y][i]
                for a in 1:lt_block_h[y][i][j], q in a:lt_block_h[y][i][j]
                    @inbounds supp_U_h[y][i][:,t_iter]=v[y][:,block_h[y][i][j][a]]+v[y][:,block_h[y][i][j][q]]
                    t_iter+=1
                end
            end
            supp_U_h[y][i]=sortslices(supp_U_h[y][i],dims=2)
            supp_U_h[y][i]=unique(supp_U_h[y][i],dims=2)
            lsupp_U_h[y][i]=size(supp_U_h[y][i],2)
            
            H[y][i]=@variable(model, [1:lsupp_U_h[y][i]])

            for a in 1:lsupp_U_h[y][i]
                for z in 1:lmon_h[W[y][i]]
                      vec=spzeros(UInt64,n)
                      vec[I[y]]=supp_U_h[y][i][:,a]+supp_h[W[y][i]][I[y],z]
                      @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,vec,n)],-coe_h[W[y][i]][z]*H[y][i][a])
                end
            end
        end
            
            
    end
        
    println("  Maximal matrix size:", maxsize)

    
    @variable(model, lambda)
                    
    cons[bfind(supp_U,lsupp_U,spzeros(UInt64,n),n)]-=lambda

    for i in 1:lmon_f
        @inbounds add_to_expression!(cons[bfind(supp_U,lsupp_U,supp_f[:,i],n)],coe_f[i])
    end 
    
    @constraint(model, cons.==0)

    @objective(model, Max, lambda)
                    
    #println(" Maximal matrix size:", maxsize)                
                    
    optimize!(model)

    opt_val = value(lambda)
    println("Termination status = ", termination_status(model))
    println("Optimal value = ",opt_val)
    
    opt_sol=[Vector{Float64}([]) for t in 1:p]
 
    return opt_val,opt_sol
end
