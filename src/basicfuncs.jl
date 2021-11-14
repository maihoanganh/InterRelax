evaluate_poly(lmon,supp,coe,a,n)=sum(coe[i]*prod(a[j]^supp[j,i] for j=1:n) for i=1:lmon)
    
    


function add_orthogonal_cons(n,lmon_g,supp_g,coe_g)

    lmon_g=[lmon_g;ones(UInt64,n)]
    alpha=zeros(UInt64,n,1)

    for j in 1:n
        alpha=zeros(UInt64,n,1)
        alpha[j]=1
        supp_g=[supp_g;[alpha]]
        coe_g=[coe_g;[ones(Float64,1)]]
    end
    
    return lmon_g,supp_g,coe_g
end

gap(a,b)=abs(a-b)/abs(b)

function get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)

    lmon_g,supp_g,coe_g=add_orthogonal_cons(n,lmon_g,supp_g,coe_g)
    m+=n
    
    @polyvar x[1:n]

    f=get_poly(lmon_f, supp_f, coe_f,x,n)
    
    g=Vector{Polynomial{true,Float64}}(undef,m)
    
    for j in 1:m
        g[j]=get_poly(lmon_g[j], supp_g[j], coe_g[j],x,n)
    end
    
    h=Vector{Polynomial{true,Float64}}(undef,l)
    
    for j in 1:l
        h[j]=get_poly(lmon_h[j], supp_h[j], coe_h[j],x,n)
    end
    
    return x,[[f];g;h]
    
end


function get_poly(lmon, supp, coe,x,n)
    f=Polynomial{true,Float64}(1.0*x[1]-x[1])
    for j=1:lmon
        f+=coe[j]*prod(x[i]^supp[i,j] for i=1:n)
    end
    return f
end


function get_info(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}};sparse::Bool=false)
    n=length(x)
    m=length(g)
    l=length(h)
            
   
    
    lmon_g=Vector{UInt64}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    
    lmon_h=Vector{UInt64}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
    
    if sparse
        supp_g=Vector{SparseMatrixCSC{UInt64}}(undef,m)
        supp_h=Vector{SparseMatrixCSC{UInt64}}(undef,l)
    else
        supp_g=Vector{Matrix{UInt64}}(undef,m)
        supp_h=Vector{Matrix{UInt64}}(undef,l)
    end
        
    dg=Vector{Int64}(undef,m)
    dh=Vector{Int64}(undef,l)



    lmon_f,supp_f,coe_f=info(f,x,n;sparse=sparse)
    
    for i in 1:m
        dg[i]=maxdegree(g[i])
        lmon_g[i],supp_g[i],coe_g[i]=info(g[i],x,n;sparse=sparse)
        
    end
                    
    for i in 1:l
        dh[i]=maxdegree(h[i])
        lmon_h[i],supp_h[i],coe_h[i]=info(h[i],x,n;sparse=sparse)
    end
    return n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh
end

function info(f,x::Vector{PolyVar{true}},n::Int64;sparse=false)
    
    mon=monomials(f)
    coe=coefficients(f)
    lmon=length(mon)
    if sparse==false
        supp=zeros(UInt64,n,lmon)
    else
        supp=spzeros(UInt64,n,lmon)
    end
    @simd for i in 1:lmon
        @simd for j in 1:n
            @inbounds supp[j,i]=DynamicPolynomials.degree(mon[i],variable(x[j]))
        end
    end
    return lmon, supp, coe
end
       


function get_basis(n::Int64,d::Int64)
    
    lb=binomial(n+d,d)
    basis=zeros(UInt64,n,lb)
    i=UInt64(0)
    t=UInt64(1)
    while i<d+1
        if basis[n,t]==i
           if i<d
              @inbounds t+=1
              @inbounds basis[1,t]=i+1
              @inbounds i+=1
           else 
                @inbounds i+=1
           end
        else 
            j=UInt64(1)
             while basis[j,t]==0
                   @inbounds j+=1
             end
             if j==1
                @inbounds t+=1
                @inbounds basis[:,t]=basis[:,t-1]
                @inbounds basis[1,t]=basis[1,t]-1
                @inbounds basis[2,t]=basis[2,t]+1
                else t+=1
                  @inbounds basis[:,t]=basis[:,t-1]
                  @inbounds basis[1,t]=basis[j,t]-1
                  @inbounds basis[j,t]=0
                  @inbounds basis[j+1,t]=basis[j+1,t]+1
             end
        end
    end
    return basis
end

#function bfind(A::Matrix{UInt64},l::Int64,a::Vector{UInt64},n::Int64)
function bfind(A,l,a,n)
    if l==0
        return 0
    end
    low=UInt64(1)
    high=l
    while low<=high
        @inbounds mid=Int(ceil(1/2*(low+high)))
        @inbounds order=comp(A[:,mid],a,n)
        if order==0
           return mid
        elseif order<0
           @inbounds low=mid+1
        else
           @inbounds high=mid-1
        end
    end
    return 0
end

#function comp(a::Vector{UInt64},b::Vector{UInt64},n::Int64)
function comp(a,b,n)
    i=UInt64(1)
    while i<=n
          if a[i]<b[i]
             return -1
          elseif a[i]>b[i]
             return 1
          else
             @inbounds i+=1
          end
    end
    if i==n+1
       return 0
    end
end

function get_info2(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}};sparse::Bool=false)
    n=length(x)
    m=length(g)
    l=length(h)
            
   
    
    lmon_g=Vector{UInt64}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    
    lmon_h=Vector{UInt64}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
    
    if sparse
        supp_g=Vector{SparseMatrixCSC{UInt64}}(undef,m)
        supp_h=Vector{SparseMatrixCSC{UInt64}}(undef,l)
    else
        supp_g=Vector{Matrix{UInt64}}(undef,m)
        supp_h=Vector{Matrix{UInt64}}(undef,l)
    end
        
    dg=Vector{Int64}(undef,m)
    dh=Vector{Int64}(undef,l)
    
    lmon_f,supp_f,coe_f=info(f,x,n;sparse=sparse)
    
    for i in 1:m
        dg[i]=maxdegree(g[i])
        lmon_g[i],supp_g[i],coe_g[i]=info(g[i],x,n;sparse=sparse)
    end
                    
    for i in 1:l
        dh[i]=maxdegree(h[i])
        lmon_h[i],supp_h[i],coe_h[i]=info(h[i],x,n;sparse=sparse)
    end
    return n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh
end


function mulpoly(n,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h)
    lsupp=lmon_g*lmon_h
    supp=Matrix{UInt64}(undef,n,lsupp)
    t=1
    for i in 1:lmon_g, j in 1:lmon_h
        supp[:,t]=supp_g[:,i]+supp_h[:,j]
        t+=1
    end
    
    
    supp=sortslices(supp,dims=2)
    supp=unique(supp,dims=2)
    lsupp=size(supp,2)
    coe=zeros(Float64,lsupp)
    
    for i in 1:lmon_g, j in 1:lmon_h
        coe[bfind(supp,lsupp,supp_g[:,i]+supp_h[:,j],n)]+=coe_g[i]*coe_h[j]
    end
    return lsupp,supp,coe
    
end
