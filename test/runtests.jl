using VPTrees
using Base.Test

# write your own tests here
@test 1 == 1

function checkintegrity{T<:Real}(tree::VPTree{T}, node::Node{T}, indices::Vector{Int})
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

checkintegrity{T<:Real}(tree::VPTree{T}, indices::Vector{Int}) = checkintegrity(tree, tree.root, indices)
