import Distances.evaluate

type PeriodicSqEuclidean{T<:AbstractFloat} <: Metric
    bw::T
end

evaluate{T<:AbstractFloat}(metric::PeriodicSqEuclidean{T}, C::Matrix, i::Int, j::Int) = distancesq(metric.bw, C, i, j)
evaluate{T<:AbstractFloat}(metric::PeriodicSqEuclidean{T}, C::Matrix, i::Int, x::T, y::T, z::T) = distancesq(metric.bw, C, i, x, y, z)

# VPTree and Node types
immutable Node
    isleaf::Bool
    index::Int
    #i::Int
    #j::Int
    distance::Float64
    left::Node
    right::Node

    Node() = new(true, zero(Int))
    Node(k::Int) = new(true, k)
    Node(k::Int, d::Float64, l::Node, r::Node) = new(false, k, d, l, r)
    #Node(k::Int, i::Int, j::Int) = new(true, k, i, j)
    #Node(k::Int, i::Int, j::Int, d::T, l::Node{T}, r::Node{T}) = new(false, k, i, j, d, l, r)
end

immutable VPTree
    n::Int
    root::Node
    metric::Metric
    points::Any
end

function partition!{T<:Real}(metric::Metric, points, indices::AbstractVector{Int}, vp::Int, d::T, i::Int, j::Int)
    while true
        @inbounds while evaluate(metric, points, vp, indices[i]) <= d
            i == j && return j
            i += 1
        end

        @inbounds while evaluate(metric, points, vp, indices[j]) > d
            i == j && return j - 1
            j -= 1
        end

        @inbounds indices[i], indices[j] = indices[j], indices[i]
    end
end

function quickselect!(metric::Metric, points, indices::AbstractVector{Int}, vp::Int, i::Int, j::Int, k::Int)
    while true
        # Choose a random element as the pivot
        r = rand(i:j)
        @inbounds d = evaluate(metric, points, vp, indices[r])
        n = partition!(metric, points, indices, vp, d, i, j)
        if n > k j = n elseif n < k i = n + 1 else return d end
    end
end

function Node(metric::Metric, points, indices::AbstractVector{Int}, i::Int, j::Int)
    j <  i && return Node()
    @inbounds j == i && return Node(indices[i])
    #@inbounds j == i && return Node{T}(indices[i], i, j)

    # Choose a random element as the vantage point
    r = rand(i:j)
    @inbounds vp = indices[r]
    @inbounds indices[i], indices[r] = indices[r], indices[i]
    i += 1

    # Find the index of the median
    m = (i + j) >> 1

    d = quickselect!(metric, points, indices, vp, i, j, m)

    left = Node(metric, points, indices, i, m)
    right = Node(metric, points, indices, m + 1, j)
    Node(vp, d, left, right)
    #Node{T}(vp, i, j, d, left, right)
end

function VPTree(metric::Metric, points, indices::AbstractVector{Int})
    n = length(indices)
    root = Node(metric, points, indices, 1, n)
    VPTree(n, root, metric, points)
end

# Compute the squared distance between two coords
function distancesq{T<:AbstractFloat}(bw::T, x1::T, y1::T, z1::T, x2::T, y2::T, z2::T)
    hbw = bw / 2
    xd = x1 - x2
    yd = y1 - y2
    zd = z1 - z2
    if xd < -hbw xd += bw elseif xd > hbw xd -= bw end
    if yd < -hbw yd += bw elseif yd > hbw yd -= bw end
    if zd < -hbw zd += bw elseif zd > hbw zd -= bw end

    return xd*xd + yd*yd + zd*zd
end

distancesq{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, i::Int, j::Int) = distancesq(bw, C[1, i], C[2, i], C[3, i], C[1, j], C[2, j], C[3, j])
distancesq{T<:AbstractFloat}(bw::T, C::AbstractMatrix{T}, i::Int, x::T, y::T, z::T) = distancesq(bw, C[1, i], C[2, i], C[3, i], x, y, z)
distancesq{T<:AbstractFloat}(bw::T, C1::AbstractMatrix{T}, C2::AbstractMatrix{T}, i1::Int, i2::Int) = distancesq(bw, C1[1, i1], C1[2, i1], C1[3, i1], C2[1, i2], C2[2, i2], C2[3, i2])
