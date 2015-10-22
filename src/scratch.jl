# Tree construction using pre-computed distances
#==function partition!{T<:Real}(v::Vector{T}, indices::Vector{Int}, i::Int, j::Int, p::T)
    while true
        while v[indices[i]] <= p
            i == j && return j
            i += 1
        end

        while v[indices[j]] > p
            i == j && return j - 1
            j -= 1
        end

        indices[i], indices[j] = indices[j], indices[i]
    end
end

function quickselect!{T<:Real}(v::Vector{T}, indices::Vector{Int}, i::Int, j::Int, k::Int)
    while true
        r = indices[rand(i:j)]
        p = v[r]
        n = partition!(v, indices, i, j, p)
        if (n > k) j = n elseif (n < k) i = n + 1 else return p end
    end
end

function nodedistances!{T<:Real}(_d::Vector{T}, _i::Vector{Int}, distance::Function, i::Int, j::Int, k::Int)
    for s in i:j
        _d[_i[s]] = distance(_i[k], _i[s])::T
    end
end

function Node{T<:Real}(_d::Vector{T}, _i::Vector{Int}, distance::Function, indices::Vector{Int}, i::Int, j::Int)
    j <  i && return Node{T}()
    j == i && return Node{T}(indices[_i[i]], i, j)

    r = rand(i:j)
    ri = _i[r]
    nodedistances!(_d, _i, distance, i, j, r)

    m = (i + j) >> 1
    d = quickselect!(_d, _i, i, j, m)

    left = Node(_d, _i, distance, indices, i, m)
    right = Node(_d, _i, distance, indices, m + 1, j)
    Node{T}(indices[ri], i, j, d, left, right)
end

function VPTree{T<:Real}(_d::Vector{T}, _i::Vector{Int}, distance::Function, indices::Vector{Int}, n::Int)
    init(_i, n)
    root = Node(_d, _i, distance, indices, 1, n)
    VPTree(n, root)
end

VPTree{T<:Real}(_d::Vector{T}, _i::Vector{Int}, distance::Function, indices::Vector{Int}) = VPTree(_d, _i, item_type, indices, length(indices))

function init(_i::Vector{Int}, n::Int)
    @inbounds @simd for i in 1:n
        _i[i] = i
    end
end

workmem{T<:Real}(::Type{VPTree}, n::Int, ::Type{T}) = return (Vector{T}(n), Vector{Int}(n)) ==#



# Tree construction using on the fly distance calculations

# Requirements for partitionu!:
#    p ∉ v
#    min(v) < p < max(v)
#==function partitionu!{T<:Real}(metric::Metric{T}, points, indices::Vector{Int}, vp::Int, d::T, i::Int, j::Int)
    while true
        @inbounds while evaluate(metric, points, vp, indices[i]) < d
            i += 1
        end

        @inbounds while evaluate(metric, points, vp, indices[j]) > d
            j -= 1
        end

        j + 1 == i && return j

        @inbounds indices[i], indices[j] = indices[j], indices[i]
    end
end

# Requirements for quickselectu!:
# j > i
function quickselectu!{T<:Real}(metric::Metric{T}, points, indices::Vector{Int}, vp::Int, i::Int, j::Int, k::Int)
    while true
        local d::T

        println([round(100*evaluate(metric, points, indices[vp], indices[s]), 2) for s in eachindex(indices)]')
        println(indices')
        println("vp: $vp\ti: $i\tj: $j\tk: $k\n")

        if i + 3 > j
            @inbounds d = (evaluate(metric, points, vp, indices[i]) + evaluate(metric, points, vp, indices[j])) / 2
        else
            r1, r2, r3 = rand(i:j-1), rand(i+1:j-1), rand(i+1:j)
            r1 == r3 && (r1 = i)
            @inbounds d = (evaluate(metric, points, vp, indices[r1])
                         + evaluate(metric, points, vp, indices[r2])
                         + evaluate(metric, points, vp, indices[r3])) / 3
        end

        n = partitionu!(metric, points, indices, vp, d, i, j)
        
        println([round(100*evaluate(metric, points, indices[vp], indices[s]), 2) for s in eachindex(indices)]')
        println(indices')
        println("vp: $vp\ti: $i\tj: $j\tk: $k\tn: $n\td: $(100d)\n\n")

        if n > k j = n elseif n < k i = n + 1 else return d end
    end
end==#

# Requirements for partitionu!:
#    p ∉ v
#    min(v) < p < max(v)
function partitionu!{T<:Real}(v::Vector{T}, indices::Vector{Int}, i::Int, j::Int, p::T)
    while true
        while v[indices[i]] < p
            i += 1
        end

        while v[indices[j]] > p
            j -= 1
        end

        j + 1 == i && return j

        indices[i], indices[j] = indices[j], indices[i]
    end
end

# Requirements for quickselectu!:
# j > i
function quickselectu!{T<:Real}(v::Vector{T}, indices::Vector{Int}, i::Int, j::Int, k::Int)
    while true
        local p::T

        if i + 3 > j
            p = (v[indices[i]] + v[indices[j]]) / 2
        else
            r1, r2, r3 = rand(i:j-1), rand(i+1:j-1), rand(i+1:j)
            r1 == r3 && (r1 = i)
            p = (v[indices[r1]] + v[indices[r2]] + v[indices[r3]]) / 3
        end

        n = partitionu!(v, indices, i, j, p)
        if (n > k) j = n elseif (n < k) i = n + 1 else return p end
    end
end

function knn{T<:Real}(tree::VPTree{T}, node::Node{T}, pq::PriorityQueue{Int, T}, τ::T, p::Int)
    d = distancesq(tree.bw, tree.points, node.index, p)

    if d < τ
        if !haskey(pq, node.index)
            dequeue!(pq)
            pq[node.index] = d
        end
        τ = peek(pq).second
    end

    node.isleaf && return τ

    μ = node.distance
    #if τ + d <= μ
    if !sumsqrtgt(τ, d, μ)
        τ = knn(tree, node.left, pq, τ, p)
    else
        τ = knn(tree, node.right, pq, τ, p)
        #if d <= τ + μ
        if sumsqrtgt(τ, μ, d)
            τ = knn(tree, node.left, pq, τ, p)
        end
    end

    return τ
end

knn{T<:Real}(tree::VPTree{T}, k::Int, pq::PriorityQueue{Int, T}, p::Int) = knn(tree, tree.root, pq, Inf, p)

function knn{T<:Real}(tree::VPTree{T}, k::Int, p::Int)
    pq = PriorityQueue(collect(-1:-1:-k), fill(Inf, k), Base.Order.Reverse)
    knn(tree, k, pq, p)
    return pq
end

knn{T<:Real}(tree::VPTree{T}, k::Int, pq::PriorityQueue{Int, T}, x::T, y::T, z::T) = knn(tree, tree.root, pq, Inf, x, y, z)

function dimmaxspread{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, x::T, y::T, z::T)
    hbw = bw / 2
    nhbw = -hbw

    n = size(C, 2)

    xmin, ymin, zmin =  Inf,  Inf,  Inf
    xmax, ymax, zmax = -Inf, -Inf, -Inf

    for i in 1:n
        xd = x - C[1, i]
        yd = y - C[2, i]
        zd = z - C[3, i]
        if xd < nhbw xd += bw elseif xd > hbw xd -= bw end
        if yd < nhbw yd += bw elseif yd > hbw yd -= bw end
        if zd < nhbw zd += bw elseif zd > hbw zd -= bw end
        if xd > xmax xmax = xd elseif xd < xmin xmin = xd end
        if yd > ymax ymax = yd elseif yd < ymin ymin = yd end
        if zd > zmax zmax = zd elseif zd < zmin zmin = zd end
    end

    return (xmax - xmin, ymax - ymin, zmax - zmin)
    return indmax((xmax - xmin, ymax - ymin, zmax - zmin))
end

function furthestpoint{T<:AbstractFloat}(cs::PointSet{T}, bw::T, x::T, y::T, z::T)
    C = cs.coords

    max_i = zero(Int)
    max_r² = zero(T)

    for i in cs.first:cs.last
        r² = distanceSqPeriodic(bw, x, y, z, C[1, i], C[2, i], C[3, i])
        if r² > max_r²
            max_r² = r²
            max_i = i
        end
    end

    return max_i, max_r²
end

function twopoints{T<:AbstractFloat}(cs::PointSet{T}, bw::T, x::T, y::T, z::T)
    C = cs.coords

    i, r² = furthestpoint(cs, bw, x, y, z)
    j_l, r²_l = furthestpoint(PointSet(cs.coords, cs.indices, cs.first, i-1), bw, C[1, i], C[2, i], C[3, i])
    j_r, r²_r = furthestpoint(PointSet(cs.coords, cs.indices, i+1, cs.last), bw, C[1, i], C[2, i], C[3, i])
    j = ifelse(r²_l > r²_r, j_l, j_r)

    return i, j
end

# Rearranges the list of points in C and returns index f such that the first f points
# are closer to point i₁ and the remainder are closer to point i₂
function splitpoints!{T<:AbstractFloat}(cs::PointSet{T}, bw::T, i₁::Int, i₂::Int)
    C = cs.coords
    f = cs.first
    l = cs.last
    while f != l
        r₁ = distanceSqPeriodic(bw, C[1, i₁], C[2, i₁], C[3, i₁], C[1, f], C[2, f], C[3, f])
        r₂ = distanceSqPeriodic(bw, C[1, i₂], C[2, i₂], C[3, i₂], C[1, f], C[2, f], C[3, f])
        if r₂ > r₁
            cs.indices[l], cs.indices[f] = cs.indices[f], cs.indices[l]
            l -= 1
        else
            f += 1
        end
    end

    return f
end

function splitpoints!{T<:AbstractFloat}(cs::PointSet{T}, x::T, y::T, z::T, bw::T, r²_split::T)
    println("Split: ", r²_split)
    C = cs.coords
    f = cs.first
    l = cs.last
    while f != l
        r² = distanceSqPeriodic(bw, C[1, f], C[2, f], C[3, f], x, y, z)
        if r² > r²_split
            cs.indices[l], cs.indices[f] = cs.indices[f], cs.indices[l]
            l -= 1
        else
            f += 1
        end
    end

    return f
end

function mindistsqperiodic{T<:AbstractFloat}(cs::PointSet{T}, bw::T, x::T, y::T, z::T)
    r²max = zero(T)
    C = cs.coords
    for i in cs.first:cs.last
        r² = distanceSqPeriodic(bw, C[1, i], C[2, i], C[3, i], x, y, z)
        if r² > r²max
            r²max = r²
        end
    end

    return r²max
end

function nearestorigin{T<:AbstractFloat}(cs::PointSet{T}, bw::T)
    C = cs.coords
    i_min = zero(Int)
    r²min = Inf

    for i in cs.first:cs.last
        r² = distanceSqPeriodic(bw, C[1, i], C[2, i], C[3, i], zero(T), zero(T), zero(T))
        if r² < r²min
            r²min = r²
            i_min = i
        end
    end

    return i_min
end

volume(r) = (4/3) * π * r^3

function randpoint{T<:AbstractFloat}(cs::PointSet{T}, bw::T)
    C = cs.coords

    i1 = rand(cs.first:cs.last)
    i2 = i1
    while i2 == i1
        i2 = rand(cs.first:cs.last)
    end

    x = (C[1, i1] + C[1, i2])/2
    y = (C[2, i1] + C[2, i2])/2
    z = (C[3, i1] + C[3, i2])/2
    x, y, z = wrap(bw, x, y, z)

    return x, y, z
end

function findgoodpoint{T<:AbstractFloat}(cs::PointSet{T}, bw::T)
    C = cs.coords
    k = length(cs) >> 1 + 1
    i = rand(cs.first:cs.last)
    #rsq = @anon j -> distanceSqPeriodic(bw, C[1, i], C[2, i], C[3, i], C[1, j], C[2, j], C[3, j])
    @inbounds lessthan = @anon (j1, j2) -> j1 < j2
    #@time sort!(cs.indices, cs.first, cs.last, PartialQuickSort(k), Base.Order.Forward)
    #return cs.indices[k]
    @time @inbounds select!(cs.indices, k, lt = isless)
end

function findgoodpoint1{T<:AbstractFloat}(cs::PointSet{T}, bw::T, m::Int, n::Int)
    C = cs.coords
    x, y, z, rsq = randpoint(cs, bw)
    println("($x\t$y\t$z)")
    println("r: $(sqrt(rsq))")
    for p in 1:m
        s = 0
        for q in 1:n
            i = rand(cs.first:cs.last)
            if distanceSqPeriodic(bw, x, y, z, C[1, i], C[2, i], C[3, i]) > rsq
                s += 1
            else
                s -= 1
            end
        end
        x = s/n
        rsq *= 2^(2x/3)
        #rsq *= 1.0 + (2/3)*x - (1/9)*x^2 # + (4/81)*x^3 - (7/243)*x^4
        println("$s/$n\t\t$(sqrt(rsq))")
    end
end

function sortpoints!{T<:AbstractFloat}(cs::PointSet{T}, bw::T, i::T)
    C = cs.coords
    f = cs.first
    l = cs.last
    while f != l
        r² = distanceSqPeriodic(bw, C[1, f], C[2, f], C[3, f], C[1, i], C[2, i], C[3, i])
        if r² > r²_split
            cs.indices[l], cs.indices[f] = cs.indices[f], cs.indices[l]
            l -= 1
        else
            f += 1
        end
    end

    return f
end

function KDTree{T<:AbstractFloat}(C::AbstractMatrix{T}, indices::AbstractVector{Int}, depth::Int)
    k = size(C, 1)
    d = mod(depth - 1, k) + 1
    l = length(indices)
    i = l >>> 1 + 1

    i <= 1 && return
    select!(indices, i, by = j -> C[d, j])
    #println(round(Int, C[indices[1:i-1]]), '\t', round(Int, C[indices[i]]), '\t', round(Int, C[indices[i+1:l]]))
    KDTree(C, sub(indices, 1:i-1), depth + 1, node << 1 + 1)
    KDTree(C, sub(indices, i+1:l), depth + 1, node << 1)
end

function Octree{T<:AbstractFloat}(C::AbstractMatrix{T}, box_width::Float64, cube_width::Float64, NUM_MAX::Int64 = 8)
    bw = box_width
    hbw = bw / 2

    o = Octree(box_width, cube_width, NUM_MAX)

    n = size(C, 2)
    for i in 1:n
        x, y, z = wrap(bw, C[1, i], C[2, i], C[3, i])
        nx = floor(Int, o.cps * ((x + hbw) / bw)) + 1
        ny = floor(Int, o.cps * ((y + hbw) / bw)) + 1
        nz = floor(Int, o.cps * ((z + hbw) / bw)) + 1

        o.num[nx, ny, nz] += 1
        p = o.num[nx, ny, nz]
        p > o.NUM_MAX && error("NUM_MAX not high enough")

        o.indices[p, nx, ny, nz] = i

        o.coords[1, p, nx, ny, nz] = x
        o.coords[2, p, nx, ny, nz] = y
        o.coords[3, p, nx, ny, nz] = z
    end
    
    return o
end

function cube{T<:Integer}(o::Octree, nx::T, ny::T, nz::T)
    cps = o.cps

    if nx <= 0 nx += cps elseif nx > cps nx -= cps end
    if ny <= 0 ny += cps elseif ny > cps ny -= cps end
    if nz <= 0 nz += cps elseif nz > cps nz -= cps end

    return nx, ny, nz
end

function cube{T<:AbstractFloat}(o::Octree, x::T, y::T, z::T)
    bw = o.box_width
    hbw = bw / 2

    nx = floor(Int, o.cps * ((x + hbw) / bw)) + 1
    ny = floor(Int, o.cps * ((y + hbw) / bw)) + 1
    nz = floor(Int, o.cps * ((z + hbw) / bw)) + 1

    return cube(o, nx, ny, nz)
end

function cube{T<:AbstractFloat}(o::Octree, coord::AbstractVector{T})
    cube(o, coord[1], coord[2], coord[3])
end

#section(x::Float64, b::Float64, n::Int) = floor(Int, n * (x + b/2) / b)

function KDTree{T<:AbstractFloat}(C::AbstractMatrix{T})
    n = size(C, 2)
    indices = collect(1:n)
    KDTree(C, indices, 1, 1)
    levels = floor(Int, log2(n)) + 1
    nodes = Vector{Int}(2^levels)

    return indices, nodes
end

function median!{T<:Real}(v::Vector{T}, indices::Vector{Int}, i::Int, j::Int)
    m = (i + j) >> 1
    return quickselect!(v, indices, i, j, m)
end

function partitionpoints!{T<:AbstractFloat}(distances::Vector{T}, bw::T, C::Matrix{T}, indices::Vector{Int}, i::Int, j::Int, k::Int)
    distances!(distances, bw, C, indices, i, j, k)
    median!(distances, indices, i, j)
end

#println("$(indices')\t$(distances[indices]')\nk: $k\tki: $ki\td: $d\n")
