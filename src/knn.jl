function knn!{T}(metric::Function, neighbors::NeighborList{T}, node::Node{T}, point_index::Int, τ::T)
    node.index == 0 && return τ

    d = metric(node.index, point_index)::T
    
    if d <= τ && d > 0.0001#&& node.index != point_index
        τ = insert_neighbor!(neighbors, node.index, d)
    end

    node.isleaf && return τ

    μ = node.radius
    if τ + d <= μ
        τ = knn!(metric, neighbors, node.inside, point_index, τ)
    else
        τ = knn!(metric, neighbors, node.outside, point_index, τ)
        if d <= τ + μ
            τ = knn!(metric, neighbors, node.inside, point_index, τ)
        end
    end

    return τ
end

function knn!{T}(neighbors::NeighborList{T}, tree::VPTree, point_index::Int)
    clear!(neighbors)
    knn!(tree.metric, neighbors, tree.root, point_index, T(Inf))
    return neighbors
end

function knn!(tree::VPTree, point_index::Int, k::Int)
    T = typeof(tree.metric(1, 1))
    neighbors = NeighborList{T}(k)
    knn!(neighbors, tree, point_index)
end
