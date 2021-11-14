function extract_optimizer(Gr::Matrix{Float64},lu0::Int64,basis_sigma0::Matrix{UInt64},opt_val::Float64,n::Int64,m::Int64,l::Int64,lmon_g::Vector{UInt64},supp_g::Vector{Matrix{UInt64}},coe_g::Vector{Vector{Float64}},lmon_h::Vector{UInt64},supp_h::Vector{Matrix{UInt64}},coe_h::Vector{Vector{Float64}},lmon_f::Int64,supp_f::Matrix{UInt64},coe_f::Vector{Float64})
    #extraction of optimizers
    V=nullspace(Gr,atol=1e-5)
    rk=size(V,2)

    println("Dimension of the null space of Gram matrix = ", rk)
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre
        V= rref_with_pivots!(Matrix(V'),1e-2)
        U=Matrix(V[1]')

        w=basis_sigma0[:,V[2]]
        N=Vector{Array{Float64}}(undef,n)
        flag=UInt8(1)
 
        
        for i in 1:n
            @inbounds kk=UInt16[]
            for j=1:size(w)[2]
                @inbounds xwj=w[:,j]
                @inbounds xwj[i]+=1
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==UInt16[]
                            @inbounds kk=t
                        else
                            @inbounds kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                @inbounds flag=0
                break
            else
                @inbounds N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            for i=1:n
                @inbounds M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                @inbounds atom=Float64[]
                for j = 1:n
                    @inbounds coordinatej=L[:,i]'*N[j]*L[:,i]
                    @inbounds push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                #println("atom ",i," = ",atom)
                println("atom ",i,":")
                @inbounds flag=1
                
                @inbounds check=evaluate_poly(lmon_f,supp_f,coe_f,atom,n)-opt_val#polynomial(f)(x => atom)-opt_val
                
                #println("  check gap of lower bound  = ",check)
                if abs(check)>1e-1
                    @inbounds flag=0
                end


                for i=1:m
                    @inbounds check=evaluate_poly(lmon_g[i],supp_g[i],coe_g[i],atom,n)#polynomial(g[i])(x => atom)
                    #println("  check inequality constraint ",i," = ",check)
                    if check<-1e-1
                        @inbounds flag=0
                    end
                end
                
                for i=1:l
                    @inbounds check=evaluate_poly(lmon_h[i],supp_h[i],coe_h[i],atom,n)
                    #println("  check equality constraint ",i," = ",check)
                    if abs(check)>1e-1
                        @inbounds flag=0
                    end
                end

                if flag ==1
                    @inbounds sol=atom
                    println("####################################")
                    #println("Optimal solution: opt_sol = ",atom)
                    println("It is an approximate optimal solution!")
                    println("####################################")
                else
                    println("It is not an approximate optimal solution!")
                end

            end
        end
    end

    return sol

end






function extract_optimizer_clique(Gr,lu0,basis_sigma0,n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h)
    #extraction of optimizers
    V=nullspace(Gr,atol=1e-5)
    rk=size(V,2)

    println("Dimension of the null space of Gram matrix = ", rk)
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre
        V= rref_with_pivots!(Matrix(V'),1e-2)
        U=Matrix(V[1]')

        w=basis_sigma0[:,V[2]]
        N=Vector{Array{Float64}}(undef,n)
        flag=UInt8(1)
 
        
        for i in 1:n
            @inbounds kk=UInt16[]
            for j=1:size(w)[2]
                @inbounds xwj=w[:,j]
                @inbounds xwj[i]+=1
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==UInt16[]
                            @inbounds kk=t
                        else
                            @inbounds kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                @inbounds flag=0
                break
            else
                @inbounds N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            for i=1:n
                @inbounds M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                @inbounds atom=Float64[]
                for j = 1:n
                    @inbounds coordinatej=L[:,i]'*N[j]*L[:,i]
                    @inbounds push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                #println("atom ",i," = ",atom)
                println("atom ",i,":")
                @inbounds flag=1
               


                for i=1:m
                    @inbounds check=evaluate_poly(lmon_g[i],supp_g[i],coe_g[i],atom,n)#polynomial(g[i])(x => atom)
                    #println("  check inequality constraint ",i," = ",check)
                    if check<-1e-1
                        @inbounds flag=0
                    end
                end
                
                for i=1:l
                    @inbounds check=evaluate_poly(lmon_h[i],supp_h[i],coe_h[i],atom,n)
                    #println("  check equality constraint ",i," = ",check)
                    if abs(check)>1e-1
                        @inbounds flag=0
                    end
                end

                if flag ==1
                    @inbounds sol=atom
                    println("############################")
                    #println("Optimal solution: opt_sol = ",atom)
                    println("It satisfies the constraints!")
                    println("############################")
                else
                    println("It does not satisfy the constraints!")
                end

            end
        end
    end

    return sol

end
