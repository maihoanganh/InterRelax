function POP_NLP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,a)
    @time begin
    model = Model(with_optimizer(Ipopt.Optimizer))
        
    set_optimizer_attribute(model, "print_level", 0)
        
    @variable(model, x[1:n])

    
    set_start_value.(x, a)

    function get_func(lmon,supp,coe)
        f=0
        for j in 1:lmon
            ind=findall(r->r>0,supp[:,j])
            lind=length(ind)
            if lind==0
                f+=coe[j]
            elseif lind==1
                f+=coe[j]*x[ind[1]]^supp[ind[1],j]
            else
                f+=coe[j]*x[ind[1]]*x[ind[2]]
            end
        end
        return f
        #return sum(coe_f[i]*prod(x[j]^supp_f[j,i] for j=1:n) for i =1:lmon_f)
    end

       
        
    @objective(model, Min, get_func(lmon_f,supp_f,coe_f))
    for i in 1:m
        @constraint(model, get_func(lmon_g[i],supp_g[i],coe_g[i])>=0)
    end
    for i in 1:l
        @constraint(model, get_func(lmon_h[i],supp_h[i],coe_h[i])==0)
    end
    optimize!(model)
    println(termination_status(model))
    opt_val=objective_value(model)
    println("opt_val=",opt_val)   
                    end
    return opt_val
end
