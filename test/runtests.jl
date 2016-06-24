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

#function exact_rangesearch(metric::


create_metric(distance::Function, points) = (i::Int, j::Int) -> distance(points, i, j)

function euclidean_sq(points::Matrix, i::Int, j::Int)
    x = points[1, i] - points[1, j]
    y = points[2, i] - points[2, j]
    z = points[3, i] - points[3, j]

    x^2 + y^2 + z^2
end

euclidean(points::Matrix, i::Int, j::Int) = sqrt(euclidean_sq(points, i, j))

function test_rangesearch(distance::Function, n::Int)
    points = randn(3, n)
    indices = collect(1:n)
    metric = create_metric(distance, points)

    @time tree = VPTree(metric, indices)
    @time rangesearch(tree, 1, 1.0)
end

test_rangesearch(euclidean, 100000)
test_rangesearch(euclidean, 100000)
