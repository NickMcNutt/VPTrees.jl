immutable Node{T}
    isleaf::Bool
    index::Int
    radius::T
    inside::Node{T}
    outside::Node{T}

    Node() = new(true, zero(Int))
    Node(k::Int) = new(true, k)
    Node(k::Int, d::T, l::Node{T}, r::Node{T}) = new(false, k, d, l, r)

    function Node(metric::Function, points, indices::AbstractVector{Int}, i::Int, j::Int)
        j <  i && return Node{T}()
        @inbounds j == i && return Node{T}(indices[i])

        # Choose a random element as the vantage point
        r = rand(i:j)
        @inbounds vp = indices[r]
        @inbounds indices[i], indices[r] = indices[r], indices[i]
        i += 1

        # Find the index of the median
        m = (i + j) >> 1

        d = quickselect!(T, metric, points, indices, vp, i, j, m)

        inside = Node{T}(metric, points, indices, i, m)
        outside = Node{T}(metric, points, indices, m + 1, j)
        Node{T}(vp, d, inside, outside)
    end
end

immutable VPTree
    num_points::Int
    root::Node
    metric::Function
    points::Any
end

function partition!(metric::Function, points, indices::AbstractVector{Int}, vp::Int, d, i::Int, j::Int)
    while true
        @inbounds while metric(points, vp, indices[i]) <= d
            i == j && return j
            i += 1
        end

        @inbounds while metric(points, vp, indices[j]) > d
            i == j && return j - 1
            j -= 1
        end

        @inbounds indices[i], indices[j] = indices[j], indices[i]
    end
end

function quickselect!{T}(::Type{T}, metric::Function, points, indices::AbstractVector{Int}, vp::Int, i::Int, j::Int, k::Int)
    while true
        # Choose a random element as the pivot
        r = rand(i:j)
        @inbounds d = metric(points, vp, indices[r])::T
        n = partition!(metric, points, indices, vp, d, i, j)

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
    VPTree(metric, points, indices)

Construct a VPTree by providing a function metric(points, i1, i2) that computes
the distance between two points in collection "points" given indices "i1" and "i2"
"""
function VPTree(metric::Function, points, indices::AbstractVector{Int})
    n = length(indices)
    T = typeof(metric(points, first(indices), first(indices)))
    root = Node{T}(metric, points, indices, 1, n)
    VPTree(n, root, metric, points)
end

function distsq(points::Matrix, i1::Int, i2::Int)
    x = points[1, i1] - points[1, i2]
    y = points[2, i1] - points[2, i2]
    z = points[3, i1] - points[3, i2]

    x^2 + y^2 + z^2
end
