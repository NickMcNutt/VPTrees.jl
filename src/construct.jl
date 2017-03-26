immutable Node{T}
    isleaf::Bool
    index::Int
    radius::T
    inside::Node{T}
    outside::Node{T}

    Node() = new(true, zero(Int))
    Node(k::Int) = new(true, k)
    Node(k::Int, d::T, l::Node{T}, r::Node{T}) = new(false, k, d, l, r)

    function Node(metric::Function, indices::Vector{Int}, i::Int, j::Int)
        j <  i && return Node{T}()
        @inbounds j == i && return Node{T}(indices[i])

        # Choose a random element as the vantage point
        r = rand(i:j)
        @inbounds vp = indices[r]
        @inbounds indices[i], indices[r] = indices[r], indices[i]
        i += 1

        # Find the index of the median
        m = (i + j) >> 1

        d = quickselect!(T, metric, indices, vp, i, j, m)::T

        inside = Node{T}(metric, indices, i, m)
        outside = Node{T}(metric, indices, m + 1, j)
        Node{T}(vp, d, inside, outside)
    end
end

immutable VPTree
    num_points::Int
    root::Node
    metric::Function
end

function partition!(metric::Function, indices::Vector{Int}, vp::Int, d, i::Int, j::Int)
    while true
        @inbounds while metric(vp, indices[i]) <= d
            i == j && return j
            i += 1
        end

        @inbounds while metric(vp, indices[j]) > d
            i == j && return j - 1
            j -= 1
        end

        @inbounds indices[i], indices[j] = indices[j], indices[i]
    end
end

function quickselect!{T}(::Type{T}, metric::Function, indices::Vector{Int}, vp::Int, i::Int, j::Int, k::Int)
    while true
        # Choose a random element as the pivot
        r = rand(i:j)
        @inbounds d = metric(vp, indices[r])::T
        n = partition!(metric, indices, vp, d, i, j)

        if n > k
            j = n
        elseif n < k
            i = n + 1
        else
            return d
        end
    end
end

"""
    VPTree(metric, indices)

Construct a VPTree by providing

1) a function metric(i::Int, j::Int) that computes the distance between point i and point j
2) a list of point indices for which the VP Tree is constructed
"""
function VPTree(metric::Function, indices::Vector{Int})
    n = length(indices)
    T = typeof(metric(first(indices), first(indices)))
    root = Node{T}(metric, indices, 1, n)
    VPTree(n, root, metric)
end

function distsq(points::Matrix, i1::Int, i2::Int)
    x = points[1, i1] - points[1, i2]
    y = points[2, i1] - points[2, i2]
    z = points[3, i1] - points[3, i2]

    x^2 + y^2 + z^2
end
