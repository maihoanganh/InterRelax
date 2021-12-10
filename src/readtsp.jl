#From https://github.com/matago/TSPLIB.jl
struct TSP
    name::AbstractString
    dimension::Integer
    weight_type::AbstractString
    weights::Matrix
    nodes::Matrix
    Dnodes::Bool
    ffx::Function
    pfx::Function
    optimal::Float64
end

const tsp_keys = ["NAME",
                "TYPE",
                "COMMENT",
                "DIMENSION",
                "EDGE_WEIGHT_TYPE",
                "EDGE_WEIGHT_FORMAT",
                "EDGE_DATA_FORMAT",
                "NODE_COORD_TYPE",
                "DISPLAY_DATA_TYPE",
                "NODE_COORD_SECTION",
                "DEPOT_SECTION",
                "DEMAND_SECTION",
                "EDGE_DATA_SECTION",
                "FIXED_EDGES_SECTION",
                "DISPLAY_DATA_SECTION",
                "TOUR_SECTION",
                "EDGE_WEIGHT_SECTION",
                "EOF"]


function readTSP(path::AbstractString)
  raw = read(path, String)
  checkEOF(raw)
  return _generateTSP(raw)
end

readTSPLIB(instance::Symbol) = readTSP(joinpath(TSPLIB95_path,string(instance)*".tsp"))

function _generateTSP(raw::AbstractString)
  _dict = keyextract(raw, tsp_keys)
  name = _dict["NAME"]
  dimension = parse(Int,_dict["DIMENSION"])
  weight_type = _dict["EDGE_WEIGHT_TYPE"]
  dxp = false

  if weight_type == "EXPLICIT" && haskey(_dict,"EDGE_WEIGHT_SECTION")
    explicits = parse.(Float64, split(_dict["EDGE_WEIGHT_SECTION"]))
    weights = explicit_weights(_dict["EDGE_WEIGHT_FORMAT"],explicits)
    #Push display data to nodes if possible
    if haskey(_dict,"DISPLAY_DATA_SECTION")
      coords = parse.(Float64, split(_dict["DISPLAY_DATA_SECTION"]))
      n_r = convert(Integer,length(coords)/dimension)
      nodes = reshape(coords,(n_r,dimension))'[:,2:end]
      dxp = true
    else
      nodes = zeros(dimension,2)
    end
  elseif haskey(_dict,"NODE_COORD_SECTION")
    coords = parse.(Float64, split(_dict["NODE_COORD_SECTION"]))
    n_r = convert(Integer,length(coords)/dimension)
    nodes = reshape(coords,(n_r,dimension))'[:,2:end]
    weights = calc_weights(_dict["EDGE_WEIGHT_TYPE"], nodes)
  end

  fFX = fullFit(weights)
  pFX = partFit(weights)

  optimal = -1
  

  TSP(name,dimension,weight_type,weights,nodes,dxp,fFX,pFX,optimal)
end

function keyextract(raw::T,ks::Array{T}) where T<:AbstractString
  pq = PriorityQueue{T,Tuple{Integer,Integer}}()
  vals = Dict{T,T}()
  for k in ks
    idx = findfirst(k,raw)
    idx != nothing && enqueue!(pq,k,extrema(idx))
  end
  while length(pq) > 1
    s_key, s_pts = peek(pq)
    dequeue!(pq)
    f_key, f_pts = peek(pq)
    rng = (s_pts[2]+1):(f_pts[1]-1)
    vals[s_key] = strip(replace(raw[rng],":"=>""))
  end
  return vals
end


function explicit_weights(key::AbstractString,data::Vector{Float64})
  w = @match key begin
    "UPPER_DIAG_ROW"  => vec2UDTbyRow(data)
    "LOWER_DIAG_ROW"  => vec2LDTbyRow(data)
    "UPPER_DIAG_COL"  => vec2UDTbyCol(data)
    "LOWER_DIAG_COL"  => vec2LDTbyCol(data)
    "UPPER_ROW"       => vec2UTbyRow(data)
    "LOWER_ROW"       => vec2LTbyRow(data)        
    "FULL_MATRIX"     => vec2FMbyRow(data)
  end
  if !in(key,["FULL_MATRIX"])
    w.+=w'
  end
  return w
end

function calc_weights(key::AbstractString, data::Matrix)
  w = @match key begin
    "EUC_2D"  => euclidian(data[:,1], data[:,2])
    "MAN_2D"  => manhattan(data[:,1], data[:,2])
    "MAX_2D"  => max_norm(data[:,1], data[:,2])
    "GEO"     => geo(data[:,1], data[:,2])
    "ATT"     => att_euclidian(data[:,1], data[:,2])
    "CEIL_2D" => ceil_euclidian(data[:,1], data[:,2])
    _         => error("Distance function type $key is not supported.")
  end

  return w
end

function checkEOF(raw::AbstractString)
  n = findlast("EOF",raw)
  if n == nothing
    throw("EOF not found")
  end
  return
end









function vec2LDTbyRow(v::AbstractVector{T}, z::T=zero(T)) where T
    n = length(v)
    s = round(Integer,(sqrt(8n+1)-1)/2)
    s*(s+1)/2 == n || error("vec2LTbyRow: length of vector is not triangular")
    k=0
    [i<=j ? (k+=1; v[k]) : z for i=1:s, j=1:s]'
end

function vec2UDTbyRow(v::AbstractVector{T}, z::T=zero(T)) where T
    n = length(v)
    s = round(Integer,(sqrt(8n+1)-1)/2)
    s*(s+1)/2 == n || error("vec2UTbyRow: length of vector is not triangular")
    k=0
    [i>=j ? (k+=1; v[k]) : z for i=1:s, j=1:s]'
end

function vec2LDTbyCol(v::AbstractVector{T}, z::T=zero(T)) where T
    n = length(v)
    s = round(Integer,(sqrt(8n+1)-1)/2)
    s*(s+1)/2 == n || error("vec2LTbyCol: length of vector is not triangular")
    k=0
    [i>=j ? (k+=1; v[k]) : z for i=1:s, j=1:s]
end

function vec2UDTbyCol(v::AbstractVector{T}, z::T=zero(T)) where T
    n = length(v)
    s = round(Integer,(sqrt(8n+1)-1)/2)
    s*(s+1)/2 == n || error("vec2UTbyCol: length of vector is not triangular")
    k=0
    [i<=j ? (k+=1; v[k]) : z for i=1:s, j=1:s]
end

function vec2UTbyRow(v::AbstractVector{T}, z::T=zero(T)) where T
    n = length(v)
    s = round(Integer,((sqrt(8n+1)-1)/2)+1)
    (s*(s+1)/2)-s == n || error("vec2UTbyRow: length of vector is not triangular")
    k=0
    [i>j ? (k+=1; v[k]) : z for i=1:s, j=1:s]'
end

function vec2LTbyRow(v::AbstractVector{T}, z::T=zero(T)) where T
    n = length(v)
    s = round(Integer,((sqrt(8n+1)-1)/2)+1)
    (s*(s+1)/2)-s == n || error("vec2LTbyRow: length of vector is not triangular")
    k=0
    [i<j ? (k+=1; v[k]) : z for i=1:s, j=1:s]'
end

function vec2FMbyRow(v::AbstractVector{T}, z::T=zero(T)) where T
    n = length(v)
    s = round(Int,sqrt(n))
    s^2 == n || error("vec2FMbyRow: length of vector is not square")
    k=0
    [(k+=1; v[k]) for i=1:s, j=1:s]
end

function findTSP(path::AbstractString)
  if isdir(path)
    syms = [Symbol(split(file,".")[1]) for file in readdir(path) if (split(file,".")[end] == "tsp")]
  else
    error("Not a valid directory")
  end
  return syms
end
                
                
#=Generator function for TSP that takes the weight Matrix
and returns a function that evaluates the fitness of a single path=#

function fullFit(costs::AbstractMatrix{Float64})
  N = size(costs,1)
  function fFit(tour::Vector{T}) where T<:Integer
    @assert length(tour) == N "Tour must be of length $N"
    @assert isperm(tour) "Not a valid tour, not a permutation"
    #distance = weights[from,to] (from,to) in tour
    distance = costs[tour[N],tour[1]]
    for i in 1:N-1
      @inbounds distance += costs[tour[i],tour[i+1]]
    end
    return distance
  end
  return fFit
end

function partFit(costs::AbstractMatrix{Float64})
  N = size(costs,1)
  function pFit(tour::Vector{T}) where T<:Integer
    n = length(tour)
    #distance = weights[from,to] (from,to) in tour
    distance = n == N ? costs[tour[N],tour[1]] : zero(Float64)
    for i in 1:n-1
      @inbounds distance += costs[tour[i],tour[i+1]]
    end
    return distance
  end
  return pFit
end

                
function euclidian(x::Vector{T},y::Vector{T}) where T<:Real
  nsz = length(x)
  dist = zeros(T,nsz,nsz)

  for i in 1:nsz, j in 1:nsz
    if i<=j
      xd = (x[i]-x[j])^2
      yd = (y[i]-y[j])^2
      dist[i,j] = dist[j,i] = round(sqrt(xd+yd),RoundNearestTiesUp)
    end
  end
  return dist
end

function manhattan(x::Vector{T},y::Vector{T}) where T<:Real
  nsz = length(x)
  dist = zeros(T, nsz, nsz)

  for i in 1:nsz, j in 1:nsz
    if i<=j
      xd = abs(x[i] - x[j])
      yd = abs(y[i] - y[j])
      dist[i, j] = dist[j, i] = round(xd + yd, RoundNearestTiesUp)
    end
  end
  return dist
end

function max_norm(x::Vector{T},y::Vector{T}) where T<:Real
  nsz = length(x)
  dist = zeros(T, nsz, nsz)

  for i in 1:nsz, j in 1:nsz
    if i<=j
      xd = abs(x[i] - x[j])
      yd = abs(y[i] - y[j])
      dist[i, j] = dist[j, i] = max(round(xd, RoundNearestTiesUp), round(yd, RoundNearestTiesUp))
    end
  end
  return dist
end

function att_euclidian(x::Vector{T},y::Vector{T}) where T<:Real
  nsz = length(x)
  dist = zeros(T,nsz,nsz)

  for i in 1:nsz, j in 1:nsz
    if i<=j
      xd = (x[i]-x[j])^2
      yd = (y[i]-y[j])^2
      dist[i,j] = dist[j,i] = ceil(sqrt((xd+yd)/10.0))
    end
  end
  return dist
end

function ceil_euclidian(x::Vector{T},y::Vector{T}) where T<:Real
  nsz = length(x)
  dist = zeros(T,nsz,nsz)

  for i in 1:nsz, j in 1:nsz
    if i<=j
      xd = (x[i]-x[j])^2
      yd = (y[i]-y[j])^2
      dist[i,j] = dist[j,i] = ceil(sqrt(xd+yd))
    end
  end
  return dist
end

function geo(x::Vector{T},y::Vector{T}) where T<:Real
  PI = 3.141592
  RRR = 6378.388
  nsz = length(x)
  dist = zeros(T,nsz,nsz)
  degs = trunc.(hcat(x,y))
  mins = hcat(x,y).-degs
  coords = PI.*(degs.+(5.0.*(mins./3.0)))./180.0
  lat = coords[:,1]
  lon = coords[:,2]

  for i in 1:nsz, j in 1:nsz
    if i<=j
      q1 = cos(lon[i]-lon[j])
      q2 = cos(lat[i]-lat[j])
      q3 = cos(lat[i]+lat[j])
      dij = RRR.*acos(0.5.*((1.0.+q1).*q2.-(1.0.-q1).*q3)).+1.0
      dist[i,j] = dist[j,i] = floor(dij)
    end
  end
  return dist
end