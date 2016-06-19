using VPTrees
using Base.Test

# write your own tests here
@test 1 == 1

#==function checkintegrity{T<:Real}(tree::VPTree{T}, node::Node{T}, indices::Vector{Int})
    i = node.i
    j = node.j

    j == i && return true

    m = (i + j) >> 1
    for k in i:j
        d = evaluate(tree.metric, tree.points, indices[k], node.index)
        ((k <= m) $ (d <= node.distance)) && return false
    end

    checkintegrity(tree, node.left, indices) && checkintegrity(tree, node.right, indices)
end

checkintegrity{T<:Real}(tree::VPTree{T}, indices::Vector{Int}) = checkintegrity(tree, tree.root, indices)==#

function distance_sq_ii(points::Matrix, i1::Int, i2::Int)
    @inbounds x = points[1, i1] - points[1, i2]
    @inbounds y = points[2, i1] - points[2, i2]
    @inbounds z = points[3, i1] - points[3, i2]

    x^2 + y^2 + z^2
end

function distance_sq_ip{T}(points::Matrix{T}, i::Int, point::NTuple{3, T})
    @inbounds x = points[1, i] - point[1]
    @inbounds y = points[2, i] - point[2]
    @inbounds z = points[3, i] - point[3]

    x^2 + y^2 + z^2
end

function inrange_test(n::Int)
    points = 10*rand(3, n) - 5
    ind = collect(1:n)

    @time tree = VPTree(distance_sq_ii, points, ind)
    neighbors = NeighborList{Float64}(n)
    @time inrange!(distance_sq_ip, neighbors, tree, 1.0, (0.0, 0.0, 0.0))

    @time for i in 1:n
        inrange!(distance_sq_ip, neighbors, tree, 1.0, (points[1, i], points[2, i], points[3, i]))
    end

    neighbors
end

inrange_test(100000)
inrange_test(100000)
