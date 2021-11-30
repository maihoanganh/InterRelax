function test()
    @polyvar x[1:2] # nonnegative variables

    f=x[1]^2+0.5*x[1]*x[2]-0.25*x[2]^2+0.75*x[1]-0.3*x[2] # the objective polynomial to minimize

    g=[1.0-sum(x.^2)] # the inequality constraints
    h=[(x[1]-1.0)*x[2]] # the equality constraints
    
    println("***Problem setting***")
    println("Number of variables: n=",2)
    m=length(g)
    println("Number of inequality constraints: m=",m)
    l=length(h)
    println("Number of equality constraints: l=",l)

    
    println()
    println("-------------------------------")
    println()
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=InterRelax.get_info(x,f,g,h,sparse=false);
    k=2
    
    TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
    println()
    println("-------------------------------")
    println()
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=InterRelax.get_info(x,f,g,h,sparse=false);
    k=1 # relaxation order
    s=3 # sparsity order
    @time RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek")
    println()
    println("-------------------------------")
    println()
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=InterRelax.get_info(x,f,g,h,sparse=false);
    k=2
    @time RelaxDense_without_multiplier(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek")
    println()
    println("-------------------------------")
    println()
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=InterRelax.get_info(x,f,g,h,sparse=true)
    k=2
    TSSOS_CS(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
    println()
    println("-------------------------------")
    println()
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=InterRelax.get_info(x,f,g,h,sparse=true)
    k=1 # relaxation order
    s=3 # sparsity order
    d=2 
    @time RelaxSparse(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,s,d,assign="min",alg="MD",minimize=true,solver="Mosek",comp_opt_sol=false)
    println()
    println("-------------------------------")
    println()
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=InterRelax.get_info(x,f,g,h,sparse=true)
    k=2
    @time RelaxSparse_without_multiplier(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,s,k,assign="min",alg="MD",minimize=true,solver="Mosek",comp_opt_sol=false);
    
end






function test_AMGM()
    
    println("***Problem setting***")

    n=3

    println("Number of variables: n=",n)

    @polyvar x[1:n]# variables

    f=sum(x)+0.0

    g=[prod(x)-1.0;3.0-sum(x)]

    m=length(g)
    println("Number of inequality constraints: m=",m)

    l=0

    h=Vector{Polynomial{true,Float64}}(undef,l)

    l=length(h)
    println("Number of equality constraints: l=",l)
    println()
    println("-------------------------------")
    println()
    
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false);
    
    
    for k in 0:20
        for s in 1:8

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false);

            @time RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek")
            println()
            println("-------------------------------")
            println()
        end
    end
end




function test_dense_POP_arbcons(data)
    
    Pb=0
    Id=0
    k=0
    s=1
    comp_opt_sol=true
    include(data*"/densePOPvar$(20)nineq$(1)neq$(0).jl")
    
    
    for n1 in [20;30]
        
        for (m1,l1) in [(1,0);(ceil(Int64, n1/5),0);(ceil(Int64, n1/5),ceil(Int64, n1/5))]
            
            Pb+=1
            
            println("Ordinal number of POP: Pb=",Pb)

            println("***Problem setting***")
            println("Number of variables: n=",n1)
            println("Number of inequality constraints: m=",m1+1)
            println("Number of equality constraints: l=",l1)
            println()
            println("-------------------------------")
            println()
            
            for k_Pu in [1;2]
                Id+=1
                println("Ordinal number of SDP: Id=",Id)
                println()
                println("-------------------------------")
                println()
                
                include(data*"/densePOPvar$(n1)nineq$(m1)neq$(l1).jl")
                k=k_Pu
                TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
                


                println()
                println("-------------------------------")
                println()

                
                include(data*"/densePOPvar$(n1)nineq$(m1)neq$(l1).jl")

                if Id in [1;2;3;5;7;8;9;10]
                    k=0
                else
                    k=1
                end
                
                if Id in [1;3;5;7;9;11]
                    s=1
                elseif Id==2
                    s=17
                elseif Id in [4;8]
                    s=20
                elseif Id==6
                    s=7
                else
                    s=31
                end
                
                if Id in [2;4;8]
                    comp_opt_sol=true
                else
                    comp_opt_sol=false
                end
                

                @time RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                    lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek",comp_opt_sol=comp_opt_sol)

                println()
                println("-------------------------------")
                println()

                
                include(data*"/densePOPvar$(n1)nineq$(m1)neq$(l1).jl")
                
                if Id in [1;2;3;5;7;8;9]
                    k=2
                else
                    k=3
                end
                
                if Id in [1;3;5;7;9;11]
                    s=1
                elseif Id==2
                    s=17
                elseif Id in [4;8]
                    s=20
                elseif Id==6
                    s=22
                elseif Id==10
                    s=32
                else
                    s=31
                end


                @time RelaxDense_without_multiplier(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek")

                println()
                println("-------------------------------")
                println("-------------------------------")
                println()
            end
        end
    end
    
end





function test_CS_POP_arbcons(data)
    n=1000
    Pb=0
    Id=0
    d=Int64(0)
    k=1
    s=1
    include(data*"/sparsePOPcliq$(10)nineq$(ceil(Int64, n/5))neq$(0).jl")
    
    for u in [10;20]
        
        for (m1,l1) in [(ceil(Int64, n/5),0);(ceil(Int64, n/5),ceil(Int64, n/5))]
            
            Pb+=1
            
            println("Ordinal number of POP: Pb=",Pb)

            println("***Problem setting***")
            println("Number of variables: n=",n)
            println("Number of inequality constraints: m=",m1+1)
            println("Number of equality constraints: l=",l1)
            println("Clique parameter: u=",u)
            println()
            println("-------------------------------")
            println()
            
            for k_Pu in [1;2]
                Id+=1
                println("Ordinal number of SDP: Id=",Id)
                println()
                println("-------------------------------")
                println()
                
                
                include(data*"/sparsePOPcliq$(u)nineq$(m1)neq$(l1).jl")
                k=k_Pu
                println("Maximal matrix size: ",binomial(u+1+k,k))
                TSSOS_CS(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)

                println()
                println("-------------------------------")
                println()

                
                include(data*"/sparsePOPcliq$(u)nineq$(m1)neq$(l1).jl")

                if Id in [1;2;3;5;6;7]
                    k=0
                else
                    k=1
                end
                
                if Id in [1;3;5]
                    s=1
                elseif Id==2
                    s=10
                elseif Id==4
                    s=12
                elseif Id==6
                    s=15
                elseif Id==7
                    s=2
                else
                    s=22
                end
                
                
                d=Int64(maximum([sum(supp_f[:,i]) for i in 1:lmon_f]))
                
                
                

                @time RelaxSparse(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,s,d,assign="min",alg="MD",minimize=true,solver="Mosek",comp_opt_sol=false)

                println()
                println("-------------------------------")
                println()

                
                include(data*"/sparsePOPcliq$(u)nineq$(m1)neq$(l1).jl")
                
                if Id in [1;2;3;5;6;7]
                    k=2
                else
                    k=3
                end
                
                if Id in [1;3;5]
                    s=1
                elseif Id==2
                    s=10
                elseif Id in [4;6]
                    s=15
                elseif Id==7
                    s=2
                else
                    s=22
                end


                


                @time RelaxSparse_without_multiplier(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,s,k,assign="min",alg="MD",minimize=true,solver="Mosek",comp_opt_sol=false);

                println()
                println("-------------------------------")
                println("-------------------------------")
                println()
            end
        end
    end
    
end






function test_PMSV(data)
    
    Pb=0
    Id=0
    k=0
    s=1
    n=16; include(data*"/mat_size$(n).jl")
    
    @polyvar x[1:n]
    f=sum(-A[i,j]*x[i]*x[j] for i=1:n for j=1:n)
    g=Vector{Polynomial{true,Float64}}([]);m=length(g)
    h=[1.0-sum(x.^2)];l=length(h)
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
    opt_val=1.0
                
    
    
    for m in [4;5;6;7]
            
        Pb+=1
        
        println("Ordinal number of POP: Pb=",Pb)
        n=m^2        
        include(data*"/mat_size$(n).jl")
        @polyvar x[1:n]
        f=sum(-A[i,j]*x[i]*x[j] for i=1:n for j=1:n)
        g=Vector{Polynomial{true,Float64}}([]);m=length(g)
        h=[1.0-sum(x.^2)];l=length(h)


        println("***Problem setting***")
        println("Number of variables: n=",n)
        println("Number of inequality constraints: m=",m+1)
        println("Number of equality constraints: l=",l)
        println()
        println("-------------------------------")
        println()

        for k_Pu in [1;2]
            Id+=1
            println("Ordinal number of SDP: Id=",Id)
            println()
            println("-------------------------------")
            println()

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=k_Pu
            opt_val=TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
                                    
            println("Approximate positive maximal singular value: sigma=",-opt_val)


            println()
            println("-------------------------------")
            println()


            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)

            k=0

            if Id==1
                s=14
            elseif Id==2
                s=17
            elseif Id==3
                s=23
            elseif Id==4
                s=26
            elseif Id==5
                s=35
            elseif Id==6
                s=37
            elseif Id==7
                s=49
            else
                s=50
            end



            @time opt_val,_=RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek",comp_opt_sol=false)
            println("Approximate positive maximal singular value: sigma=",-opt_val)


            println()
            println("-------------------------------")
            println("-------------------------------")
            println()
        end
    end
    
end
                        
                        
                        

function test_compute_stability_number_of_graph_random(data)
    
    Pb=0
    Id=0
    k=0
    s=1
    n=10; include(data*"/mat_stability_size$(n).jl")
    
    @polyvar x[1:n]
    f=sum(A[i,j]*x[i]*x[j] for i=1:n for j=1:n)
    g=Vector{Polynomial{true,Float64}}([]);m=length(g)
    h=[1.0-sum(x)];l=length(h)
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
    opt_val=1.0
                
    
    
    for n1 in [10;15;20;25]
            
        Pb+=1
        
        println("Ordinal number of POP: Pb=",Pb)    
        n=n1
        include(data*"/mat_stability_size$(n).jl")
        @polyvar x[1:n]
        f=sum(A[i,j]*x[i]*x[j] for i=1:n for j=1:n)
        g=Vector{Polynomial{true,Float64}}([]);m=length(g)
        h=[1.0-sum(x)];l=length(h)
                                                       


        println("***Problem setting***")
        println("Number of variables: n=",n)
        println("Number of inequality constraints: m=",m+1)
        println("Number of equality constraints: l=",l)
        println()
        println("-------------------------------")
        println()

        for k_Pu in [1;2]
            Id+=1
            println("Ordinal number of SDP: Id=",Id)
            println()
            println("-------------------------------")
            println()

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=k_Pu
            opt_val=TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
            println("Approximate stability number: alpha=",1/opt_val)


            println()
            println("-------------------------------")
            println()


            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            if Id in [1;2;3;4;5;7]
                k=0
            else
                k=1
            end
                                                            
            if Id==1
                s=6
            elseif Id==2
                s=11
            elseif Id==3
                s=10
            elseif Id==4
                s=16
            elseif Id==5
                s=21
            elseif Id==6
                s=21
            elseif Id==7
                s=26
            else
                s=26
            end



            @time opt_val,_=RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek",comp_opt_sol=false)

            println("Approximate stability number: alpha=",1/opt_val)

            println()
            println("-------------------------------")
            println("-------------------------------")
            println()
        end
    end
    
end
                                                
                                                
                                                
                                                
function test_compute_stability_number_of_graph(data)
    
    Pb="GD02_a"
    Id=0
    k=0
    s=1

    A = MatrixMarket.mmread(data*"/GD02_a.mtx")

    n=size(A,1)
    for i in 1:n
        A[i,i]=1
    end
                                                    
    @polyvar x[1:n]
    f=sum(A[i,j]*x[i]*x[j] for i=1:n for j=1:n)*1.0
    g=Vector{Polynomial{true,Float64}}([]);m=length(g)
    h=[1.0-sum(x)];l=length(h)
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
    opt_val=1.0
                
    
    
    for n1 in 1:5
                                                                    
                                                                    
        if n1==1
            A = MatrixMarket.mmread(data*"/GD02_a.mtx")
            Pb="GD02_a"
        elseif n1==2
            A = MatrixMarket.mmread(data*"/johnson8-2-4.mtx")
            Pb="johnson8-2-4"
        elseif n1==3
            A = MatrixMarket.mmread(data*"/johnson8-4-4.mtx")
            Pb="johnson8-4-4"
        elseif n1==4
            A = MatrixMarket.mmread(data*"/hamming6-2.mtx")
            Pb="hamming6-2"
        else
            A = MatrixMarket.mmread(data*"/hamming6-4.mtx")
            Pb="hamming6-4"
        end
                                                                    
        n=size(A,1)
        for i in 1:n
            A[i,i]=1
        end                                                            

        
        println("Ordinal number of POP: Pb=",Pb)    
                                                                        
        @polyvar x[1:n]
        f=sum(A[i,j]*x[i]*x[j] for i=1:n for j=1:n)*1.0
        g=Vector{Polynomial{true,Float64}}([]);m=length(g)
        h=[1.0-sum(x)];l=length(h)
                                                       


        println("***Problem setting***")
        println("Number of variables: n=",n)
        println("Number of inequality constraints: m=",m+1)
        println("Number of equality constraints: l=",l)
        println()
        println("-------------------------------")
        println()

        for k_Pu in [1;2]
            Id+=1
            println("Ordinal number of SDP: Id=",Id)
            println()
            println("-------------------------------")
            println()

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=k_Pu
            opt_val=TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
            println("Approximate stability number: alpha=",1/opt_val)


            println()
            println("-------------------------------")
            println()


            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=0
                                                            
            if Id==1
                s=12
            elseif Id==2
                s=25
            elseif Id==3
                s=20
            elseif Id==4
                s=30
            elseif Id==5
                s=60
            elseif Id==6
                s=72
            elseif Id==7
                s=40
            elseif Id==8
                s=66
            elseif Id==9
                s=40
            else
                s=60
            end



            @time opt_val,_=RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek",comp_opt_sol=false)

            println("Approximate stability number: alpha=",1/opt_val)

            println()
            println("-------------------------------")
            println("-------------------------------")
            println()
        end
    end
    
end
                                                                        
                                                                        
                                                                        
function test_deciding_copositivity(data)
    
    Pb=0
    Id=0
    k=0
    s=1
    n=10; include(data*"/mat_copositivity_size$(n).jl")
    
    @polyvar x[1:n]
    f=sum(A[i,j]*x[i]*x[j] for i=1:n for j=1:n)
    g=Vector{Polynomial{true,Float64}}([]);m=length(g)
    h=[1.0-sum(x)];l=length(h)
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
    opt_val=1.0
                
    
    
    for n1 in [10;15;20;25]
            
        Pb+=1
        
        println("Ordinal number of POP: Pb=",Pb)    
        n=n1
        include(data*"/mat_copositivity_size$(n).jl")
        @polyvar x[1:n]
        f=sum(A[i,j]*x[i]*x[j] for i=1:n for j=1:n)
        g=Vector{Polynomial{true,Float64}}([]);m=length(g)
        h=[1.0-sum(x)];l=length(h)
                                                       


        println("***Problem setting***")
        println("Number of variables: n=",n)
        println("Number of inequality constraints: m=",m+1)
        println("Number of equality constraints: l=",l)
        println()
        println("-------------------------------")
        println()

        for k_Pu in [1;2]
            Id+=1
            println("Ordinal number of SDP: Id=",Id)
            println()
            println("-------------------------------")
            println()

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=k_Pu
            opt_val=TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)


            println()
            println("-------------------------------")
            println()


            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=0
                                                            
            if Id==1
                s=4
            elseif Id==2
                s=8
            elseif Id==3
                s=8
            elseif Id==4
                s=13
            elseif Id==5
                s=10
            elseif Id==6
                s=20
            elseif Id==7
                s=10
            else
                s=19
            end



            @time opt_val,_=RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek",comp_opt_sol=true)


            println()
            println("-------------------------------")
            println("-------------------------------")
            println()
        end
    end
    
end
                                                                                                
                                                                                                
                                                                                
function test_deciding_nonegativity(data)
    
    Pb=0
    Id=0
    k=0
    s=1
    n=5; include(data*"/vec_nonneg_var$(n).jl")
    
    @polyvar x[1:n]
    v=reverse(monomials(x,4))
    f=c'*v
    g=Vector{Polynomial{true,Float64}}([]);m=length(g)
    h=[1.0-sum(x)];l=length(h)
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
    opt_val=1.0
                
    
    
    for n1 in [5;10;15]
            
        Pb+=1
        
        println("Ordinal number of POP: Pb=",Pb)    
        n=n1
        include(data*"/vec_nonneg_var$(n).jl")
    
        @polyvar x[1:n]
        v=reverse(monomials(x,4))
        f=c'*v
        g=Vector{Polynomial{true,Float64}}([]);m=length(g)
        h=[1.0-sum(x)];l=length(h)



        println("***Problem setting***")
        println("Number of variables: n=",n)
        println("Number of inequality constraints: m=",m+1)
        println("Number of equality constraints: l=",l)
        println()
        println("-------------------------------")
        println()

        for k_Pu in [2;3]
            Id+=1
            println("Ordinal number of SDP: Id=",Id)
            println()
            println("-------------------------------")
            println()

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=k_Pu
            opt_val=TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)


            println()
            println("-------------------------------")
            println()


            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=0
                                                            
            if Id==1
                s=4
            elseif Id==2
                s=8
            elseif Id==3
                s=5
            elseif Id==4
                s=11
            elseif Id==5
                s=20
            else
                s=44
            end



            @time opt_val,_=RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek",comp_opt_sol=true)


            println()
            println("-------------------------------")
            println("-------------------------------")
            println()
        end
    end
    
end
                                                                                                
                                                                                                
                                                                                                
                                                                                                
function test_dense_POP_binary_constr_random(data)
    
    Pb=0
    Id=0
    k=0
    s=1
    n=10; include(data*"/vec_binary_var$(n).jl")
    
    @polyvar x[1:n]
    v=reverse(monomials(x,0:2))
    f=c'*v
    g=Vector{Polynomial{true,Float64}}([]);m=length(g)
    h=[x[j]*(1.0-x[j]) for j=1:n];l=length(h)
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
    opt_val=1.0
                
    
    
    for n1 in [10;20;30]
            
        Pb+=1
        
        println("Ordinal number of POP: Pb=",Pb)    
        n=n1
        include(data*"/vec_binary_var$(n).jl")
    
        @polyvar x[1:n]
        v=reverse(monomials(x,0:2))
        f=c'*v
        g=Vector{Polynomial{true,Float64}}([]);m=length(g)
        h=[x[j]*(1.0-x[j]) for j=1:n];l=length(h)



        println("***Problem setting***")
        println("Number of variables: n=",n)
        println("Number of inequality constraints: m=",m+1)
        println("Number of equality constraints: l=",l)
        println()
        println("-------------------------------")
        println()

        for k_Pu in [1;2]
            Id+=1
            println("Ordinal number of SDP: Id=",Id)
            println()
            println("-------------------------------")
            println()

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=k_Pu
            opt_val=TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)


            println()
            println("-------------------------------")
            println()


            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=false)
            k=1
                                                            
            if Id==1
                s=7
            elseif Id==2
                s=11
            elseif Id==3
                s=15
            elseif Id==4
                s=21
            elseif Id==5
                s=20
            else
                s=31
            end



            @time opt_val,_=RelaxDense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,
                lmon_f,supp_f,coe_f,dg,dh,k,s,solver="Mosek",comp_opt_sol=false)


            println()
            println("-------------------------------")
            println("-------------------------------")
            println()
        end
    end
    
end






function test_CertifyNN(data)
    n=1000
    Pb=0
    Id=0
    d=Int64(0)
    k=1
    s=1
    
    
    test_x, test_y = MNIST.testdata();
    vars = matread(data*"/deq_MNIST_SingleFcNet_m=20.mat");
    A = vars["A"]; B = vars["B"]; U = vars["U"]; u = vars["u"]; C = vars["C"];
    p = size(U, 1); m = 20; W = (1-m)*Matrix(I(p))-A'*A+B-B';

    U=U[:,1:p]

    me = 0.1307; std = 0.3081; eps = 0.1/std;
    k = 1; x0 = (vec(test_x[:,:,k]).-me)./std; 


    p0=size(U, 2)
    @polyvar x_bar[1:p0] z[1:p]# variables
    n=p0+p
    x=[x_bar;z]

    f=-C[1,:]'*z

    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])

    pol=1.0*z[1]

    for j in 1:p
        pol=z[j]-sum(U[j,r]*(x_bar[r]+x0[r]-eps) for r=1:p0)-u[j]
        append!(g,[pol])
        append!(h,[z[j]*pol])
    end

    for t in 1:p0
        append!(g,[2*eps-x_bar[t]])
    end

    
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=true)
    
    
    for i in [1;2;3]
        
        f=-C[i,:]'*z
        
        n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=true)

        println("Ordinal number of POP: i=",i)

        println("***Problem setting***")
        println("Number of variables: n=",n)
        println("Number of inequality constraints: m=",m+1)
        println("Number of equality constraints: l=",l)
        println()
        println("-------------------------------")
        println()

        for k_Pu in [1;2]
            Id+=1
            println("Ordinal number of SDP: Id=",Id)
            println()
            println("-------------------------------")
            println()

            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=true)
            k=k_Pu
            println("Maximal matrix size: ",binomial(88+k,k))
            TSSOS_CS(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)

            println()
            println("-------------------------------")
            println()


            n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=get_info(x,f,g,h,sparse=true)
            
            k=0

            if k_Pu==1
                s=20
            else
                s=90
            end


            d=Int64(maximum([sum(supp_f[:,i]) for i in 1:lmon_f]))+1




            @time RelaxSparse(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh,k,s,d,assign="min",alg="MD",minimize=true,solver="Mosek",comp_opt_sol=false)





            println()
            println("-------------------------------")
            println("-------------------------------")
            println()
        end
    end
    
end

