function TSSOS_Dense(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
    println("**Semidefinite relaxation based on Putinar's Positivstellensatz**")
    println("Relaxation order: k=",k)
    println("Maximal matrix size: ",binomial(n+k,k))
    
    
    
    opt=0.0
    vars,pop=get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)
    
    if binomial(n+k,k)>500
        println("Out of memory!")
    else
        try
            @time opt,sol,data=TSSOS.cs_tssos_first(pop,vars,k,numeq=l,CS=false,TS=false,solver="Mosek")
        catch
            println("Out of memory!")
        end
    end
    return opt
end


function TSSOS_CS(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,k)
    #println("**Semidefinite relaxation based on Putinar's Positivstellensatz**")
    #println("Relaxation order: k=",k)
   
    
    opt=0.0
    vars,pop=get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)
    try
        @time opt,sol,data=TSSOS.cs_tssos_first(pop,vars,k,numeq=l,CS="MD",TS=false,solver="Mosek")
    catch
        println("Out of memory!")
    end
    return opt
end
    