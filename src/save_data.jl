function save_info_sparsePOP(randx,u,n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh)
    data="/home/hoanganh/Desktop/math-topics/InterruptedRelax/codes/dataPOP2" # path of data
    output_file = open(data*"/sparsePOPcliq$(u)nineq$(m)neq$(l).jl","w")

    write(output_file, "# We are saving information of POP here. \n \n")
    
    write(output_file, "randx = ")
    show(output_file, randx)
    write(output_file, "; \n \n")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "m = ")
    show(output_file, m)
    write(output_file, "; \n \n")
    
    write(output_file, "l = ")
    show(output_file, l)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_g = ")
    show(output_file, lmon_g)
    write(output_file, "; \n \n")
    
    I_g=Vector{Vector{UInt64}}(undef,m)
    J_g=Vector{Vector{UInt64}}(undef,m)
    V_g=Vector{Vector{UInt64}}(undef,m)
    
    for i in 1:m
        I_g[i],J_g[i],V_g[i]=findnz(supp_g[i])
    end
    
    write(output_file, "I_g = ")
    show(output_file, I_g)
    write(output_file, "; \n \n")
    
    write(output_file, "J_g = ")
    show(output_file, J_g)
    write(output_file, "; \n \n")
    
    write(output_file, "V_g = ")
    show(output_file, V_g)
    write(output_file, "; \n \n")
        
    
    write(output_file, "supp_g = Vector{SparseMatrixCSC{UInt64}}([sparse(I_g[i],J_g[i],V_g[i],n,lmon_g[i]) for i in 1:m])")
    write(output_file, "; \n \n")
    
    write(output_file, "coe_g = ")
    show(output_file, coe_g)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_h = ")
    show(output_file, lmon_h)
    write(output_file, "; \n \n")
    
    I_h=Vector{Vector{UInt64}}(undef,l)
    J_h=Vector{Vector{UInt64}}(undef,l)
    V_h=Vector{Vector{UInt64}}(undef,l)
    
    for i in 1:l
        I_h[i],J_h[i],V_h[i]=findnz(supp_h[i])
    end
    
    write(output_file, "I_h = ")
    show(output_file, I_h)
    write(output_file, "; \n \n")
    
    write(output_file, "J_h = ")
    show(output_file, J_h)
    write(output_file, "; \n \n")
    
    write(output_file, "V_h = ")
    show(output_file, V_h)
    write(output_file, "; \n \n")
        
    
    write(output_file, "supp_h = Vector{SparseMatrixCSC{UInt64}}([sparse(I_h[i],J_h[i],V_h[i],n,lmon_h[i]) for i in 1:l])")
    write(output_file, "; \n \n")
    
    write(output_file, "coe_h = ")
    show(output_file, coe_h)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_f = ")
    show(output_file, lmon_f)
    write(output_file, "; \n \n")
    
    I_f,J_f,V_f=findnz(supp_f)
    
    write(output_file, "I_f = ")
    show(output_file, I_f)
    write(output_file, "; \n \n")
    
    write(output_file, "J_f = ")
    show(output_file, J_f)
    write(output_file, "; \n \n")
    
    write(output_file, "V_f = ")
    show(output_file, V_f)
    write(output_file, "; \n \n")
        
    
    write(output_file, "supp_f = sparse(I_f,J_f,V_f,n,lmon_f)")
    write(output_file, "; \n \n")
    
    write(output_file, "coe_f = ")
    show(output_file, coe_f)
    write(output_file, "; \n \n")
    
    write(output_file, "dg = ")
    show(output_file, dg)
    write(output_file, "; \n \n")
    
    write(output_file, "dh = ")
    show(output_file, dh)
    write(output_file, "; \n \n")
    
    close(output_file)
    
end



function save_info_densePOP(randx,n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh)
    data="/home/hoanganh/Desktop/math-topics/InterruptedRelax/codes/dataPOP2" # path of data
    output_file = open(data*"/densePOPvar$(n)nineq$(m)neq$(l).jl","w")

    write(output_file, "# We are saving information of POP here. \n \n")
    
    write(output_file, "randx = ")
    show(output_file, randx)
    write(output_file, "; \n \n")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "m = ")
    show(output_file, m)
    write(output_file, "; \n \n")
    
    write(output_file, "l = ")
    show(output_file, l)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_g = ")
    show(output_file, lmon_g)
    write(output_file, "; \n \n")
    
    write(output_file, "supp_g = ")
    show(output_file, supp_g)
    write(output_file, "; \n \n")
    
    write(output_file, "coe_g = ")
    show(output_file, coe_g)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_h = ")
    show(output_file, lmon_h)
    write(output_file, "; \n \n")
    
    write(output_file, "supp_h = ")
    show(output_file, supp_h)
    write(output_file, "; \n \n")
    
    write(output_file, "coe_h = ")
    show(output_file, coe_h)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_f = ")
    show(output_file, lmon_f)
    write(output_file, "; \n \n")
    
    write(output_file, "supp_f = ")
    show(output_file, supp_f)
    write(output_file, "; \n \n")
    
    write(output_file, "coe_f = ")
    show(output_file, coe_f)
    write(output_file, "; \n \n")
    
    write(output_file, "dg = ")
    show(output_file, dg)
    write(output_file, "; \n \n")
    
    write(output_file, "dh = ")
    show(output_file, dh)
    write(output_file, "; \n \n")
    
    close(output_file)
    
end






function save_info_mat(n,A)
    data="/home/hoanganh/Desktop/math-topics/InterruptedRelax/codes/dataPOP2" # path of data
    output_file = open(data*"/mat_size$(n).jl","w")

    write(output_file, "# We are saving information of matrix here. \n \n")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "A = ")
    show(output_file, A)
    write(output_file, "; \n \n")
    
  
    
    close(output_file)
    
end


function save_info_mat_stability(n,A)
    data="/home/hoanganh/Desktop/math-topics/InterruptedRelax/codes/dataPOP2" # path of data
    output_file = open(data*"/mat_stability_size$(n).jl","w")

    write(output_file, "# We are saving information of matrix here. \n \n")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "A = ")
    show(output_file, A)
    write(output_file, "; \n \n")
    
  
    
    close(output_file)
    
end

function save_info_mat_copositivity(n,A)
    data="/home/hoanganh/Desktop/math-topics/InterruptedRelax/codes/dataPOP2" # path of data
    output_file = open(data*"/mat_copositivity_size$(n).jl","w")

    write(output_file, "# We are saving information of matrix here. \n \n")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "A = ")
    show(output_file, A)
    write(output_file, "; \n \n")
    
  
    
    close(output_file)
    
end

function save_info_vec_nonneg(n,c)
    data="/home/hoanganh/Desktop/math-topics/InterruptedRelax/codes/dataPOP2" # path of data
    output_file = open(data*"/vec_nonneg_var$(n).jl","w")

    write(output_file, "# We are saving information of matrix here. \n \n")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "c = ")
    show(output_file, c)
    write(output_file, "; \n \n")
    
  
    
    close(output_file)
    
end


function save_info_vec_binary(n,c)
    data="/home/hoanganh/Desktop/math-topics/InterruptedRelax/codes/dataPOP2" # path of data
    output_file = open(data*"/vec_binary_var$(n).jl","w")

    write(output_file, "# We are saving information of matrix here. \n \n")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "c = ")
    show(output_file, c)
    write(output_file, "; \n \n")
    
  
    
    close(output_file)
    
end