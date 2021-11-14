function pop_opf(data::Dict{String, Any}; normalize = true)

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    model = PolyModel()   
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])

    vr = Dict()
    vi = Dict()

    for key in keys(ref[:bus])
        # real part voltage variables
        vr[key] = PolyPowerModels.new_polyvar("vr"*string(key))
        add_constraint!(model, ( vr[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vr[key] ), GT; normalize = normalize)

        # imaginary part voltage variables
        if haskey(ref[:ref_buses], key)
            vi[key] = 0.0
        else
            vi[key] = PolyPowerModels.new_polyvar("vi"*string(key))
            add_constraint!(model, ( vi[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vi[key] ), GT; normalize = normalize)
        end

        # voltage magnitude constraints
        add_constraint!(model, ref[:bus][key]["vmin"]^2, LT, vr[key]^2 + vi[key]^2; normalize = normalize )
        add_constraint!(model, vr[key]^2 + vi[key]^2, LT, ref[:bus][key]["vmax"]^2; normalize = normalize )

    end

    p = Dict()
    q = Dict()

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        vr_fr = vr[branch["f_bus"]]
        vr_to = vr[branch["t_bus"]]
        vi_fr = vi[branch["f_bus"]]
        vi_to = vi[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        p[f_idx] = PolyPowerModels.new_polyvar("p"*string(f_idx))
        add_constraint!( model, p[f_idx], EQ, 
                        (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) 
                        + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to)
                        + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to);
                        normalize = normalize )

        q[f_idx] = PolyPowerModels.new_polyvar("q"*string(f_idx))  
        add_constraint!( model, q[f_idx], EQ, 
                        -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) 
                        - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) 
                        + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to);
                        normalize = normalize )

        p[t_idx] = PolyPowerModels.new_polyvar("p"*string(t_idx))
        add_constraint!( model, p[t_idx], EQ, 
                        (g+g_to)*(vr_to^2 + vi_to^2) 
                        + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) 
                        + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to));
                        normalize = normalize )

        q[t_idx] = PolyPowerModels.new_polyvar("q"*string(t_idx))
        add_constraint!( model, q[t_idx], EQ, 
                        -(b+b_to)*(vr_to^2 + vi_to^2) 
                        - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to)
                        + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to));
                        normalize = normalize )

        # angle differences
        add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), LT, tan(branch["angmax"])*(vr_fr*vr_to + vi_fr*vi_to); normalize = normalize )
        add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), GT, tan(branch["angmin"])*(vr_fr*vr_to + vi_fr*vi_to); normalize = normalize )

        # thermal limits
        add_constraint!( model, p[f_idx]^2 + q[f_idx]^2, LT, branch["rate_a"]^2; normalize = normalize )
        add_constraint!( model, p[t_idx]^2 + q[t_idx]^2, LT, branch["rate_a"]^2; normalize = normalize )
    end

    for (l,i,j) in ref[:arcs]
        add_constraint!( model, -ref[:branch][l]["rate_a"], LT, p[(l,i,j)]; normalize = normalize )
        add_constraint!( model,  p[(l,i,j)], LT, ref[:branch][l]["rate_a"]; normalize = normalize )
        add_constraint!( model, -ref[:branch][l]["rate_a"], LT, q[(l,i,j)]; normalize = normalize )
        add_constraint!( model,  q[(l,i,j)], LT, ref[:branch][l]["rate_a"]; normalize = normalize )
    end


    pg = Dict()
    qg = Dict()

    for (i,bus) in ref[:bus]

        # active/reactive power

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        if isempty(ref[:bus_gens][i])
            add_constraint!(model,PolyPowerModels.fl_sum(p[a] for a in ref[:bus_arcs][i]),EQ,-PolyPowerModels.fl_sum(load["pd"] for load in bus_loads) - PolyPowerModels.fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize )
            add_constraint!(model,PolyPowerModels.fl_sum(q[a] for a in ref[:bus_arcs][i]),EQ,-PolyPowerModels.fl_sum(load["qd"] for load in bus_loads) + PolyPowerModels.fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize )
        else
            for gen_id in ref[:bus_gens][i]
                pg[gen_id] = PolyPowerModels.new_polyvar("pg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["pmin"], LT, pg[gen_id]; normalize = normalize )
                add_constraint!( model, pg[gen_id], LT, ref[:gen][gen_id]["pmax"]; normalize = normalize )

                qg[gen_id] = PolyPowerModels.new_polyvar("qg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["qmin"], LT, qg[gen_id]; normalize = normalize )
                add_constraint!( model, qg[gen_id], LT, ref[:gen][gen_id]["qmax"]; normalize = normalize )

            end
            add_constraint!( model,
                            PolyPowerModels.fl_sum(p[a] for a in ref[:bus_arcs][i]), EQ, 
                            PolyPowerModels.fl_sum(pg[g] for g in ref[:bus_gens][i]) -
                            PolyPowerModels.fl_sum(load["pd"] for load in bus_loads) -
                            PolyPowerModels.fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize
                           )  
            add_constraint!( model,
                            PolyPowerModels.fl_sum(q[a] for a in ref[:bus_arcs][i]), EQ,
                            PolyPowerModels.fl_sum(qg[g] for g in ref[:bus_gens][i]) -
                            PolyPowerModels.fl_sum(load["qd"] for load in bus_loads) +
                            PolyPowerModels.fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize
                           )
        end
    end

    # objective
    set_objective!( model, MIN, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )

    var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg)

    return PolyPowerModel(model, data, ref, var, Dict())
end


function get_POP_OPF_normal(data::Dict{String, Any})
    pm = pop_opf(data,normalize = false)
    #pm = pop_opf_deg2(data, normalize = true)
    
    f=objective_function(pm)
    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])
    
    cons=constraints(pm)
    for j in 1:length(cons)
        if sense(cons[j])==PolyPowerModels.EQ_sense()
            push!(h,constraint_function(cons[j]))
        elseif sense(cons[j])==PolyPowerModels.LT_sense()
            push!(g,-constraint_function(cons[j]))
        else
            push!(g,constraint_function(cons[j]))
        end
    end
    
    x=unique(vcat(variables.([[f];g;h])...))
    
    starpoint=ones(Float64,length(x))
    
    return x,f,g,h,starpoint
end


function get_POP_OPF(data::Dict{String, Any})
    pm = pop_opf(data,normalize = true)
    #pm = pop_opf_deg2(data, normalize = true)
    
    f=objective_function(pm)
    g=Vector{Polynomial{true,Float64}}([])
    h=Vector{Polynomial{true,Float64}}([])
    
    cons=constraints(pm)
    for j in 1:length(cons)
        if sense(cons[j])==PolyPowerModels.EQ_sense()
            push!(h,constraint_function(cons[j]))
        elseif sense(cons[j])==PolyPowerModels.LT_sense()
            push!(g,-constraint_function(cons[j]))
        else
            push!(g,constraint_function(cons[j]))
        end
    end
    
    x=unique(vcat(variables.([[f];g;h])...))
    
    #x = sort!(union(variables(f),variables.([g;h])...), rev = true)
    
    
    
    
    mon=monomials(g[1])
    coe=coefficients(g[1])
    lmon=length(mon)

    ind=1

    n=length(x)
    m=length(g)
    l=length(h)
    lb=Vector{Float64}([-Inf for j in 1:n])
    for i in 1:m
        mon=monomials(g[i])
        coe=coefficients(g[i])
        lmon=length(mon)
        if mon[lmon]==1
            for j in 1:lmon-1
                ind=findfirst(y->y==mon[j],x)
                if ind != nothing && coe[j]>0
                    lb[ind]=maximum([lb[ind];-coe[lmon]/coe[j]])
                else
                    ind=findfirst(y->y==mon[j],x.^2)
                    if ind != nothing && coe[j]<0
                        lb[ind]=maximum([lb[ind];-sqrt(-coe[lmon]/coe[j])])
                    end
                end
            end
        end
    end
    
    nlb=findall(y->y==-Inf,lb)
    @polyvar xadd[1:length(nlb)]

    y=x+lb
    y[nlb]=x[nlb]-xadd

    f=f(x=>y)
    for j in 1:m
        g[j]=g[j](x=>y)
    end
    
    for j in 1:l
        h[j]=h[j](x=>y)
    end
    
    rmind=Vector{UInt64}([])
    for i in 1:m
        if length(monomials(g[i]))==1
            append!(rmind,i)
        end
    end
    g=g[setdiff(1:m,rmind)]
    
    x=unique(vcat(variables.([[f];g;h])...))
    
    starpoint=ones(Float64,length(x))
    indzero=findall(y->y in xadd, x)
    starpoint[indzero]=zeros(Float64,length(indzero))
    
    
    
    
    
    return x,f,g,h,starpoint
end


