function rangesearch!{T}(metric::Function, neighbors::NeighborList{T}, node::Node{T}, point_index::Int, r::T)
    node.index == 0 && return

    d = metric(point_index, node.index)::T
    d <= r && add_neighbor!(neighbors, node.index, d)

    node.isleaf && return

    μ = node.radius
    if d + μ <= r
        add_neighbors!(metric, neighbors, node.inside, point_index)
    elseif d - μ <= r
        rangesearch!(metric, neighbors, node.inside, point_index, r)
    end

    if μ - d <= r
        rangesearch!(metric, neighbors, node.outside, point_index, r)
    end
end

function rangesearch_sq!{T}(metric::Function, neighbors::NeighborList{T}, node::Node{T}, point_index::Int, r::T)
    node.index == 0 && return

    d = metric(point_index, node.index)::T
    d <= r && add_neighbor!(neighbors, node.index, d)

    node.isleaf && return

    μ = node.radius
    # d + μ <= r
    if !sumsqrtgt(d, μ, r)
        add_neighbors!(metric, neighbors, node.inside, point_index)
    # d - μ <= r
    elseif sumsqrtgt(r, μ, d)
        rangesearch_sq!(metric, neighbors, node.inside, point_index, r)
    end

    # μ - d <= r
    if sumsqrtgt(r, d, μ)
        rangesearch_sq!(metric, neighbors, node.outside, point_index, r)
    end
end

function rangesearch!{T}(neighbors::NeighborList{T}, tree::VPTree, point_index::Int, r::T)
    clear!(neighbors)
    rangesearch!(tree.metric, neighbors, tree.root, point_index, r)
end

function rangesearch_sq!{T}(neighbors::NeighborList{T}, tree::VPTree, point_index::Int, r::T)
    clear!(neighbors)
    rangesearch_sq!(tree.metric, neighbors, tree.root, point_index, r)
end

function rangesearch{T}(tree::VPTree, point_index::Int, r::T)
    neighbors = NeighborList{T}(tree.num_points)
    rangesearch!(neighbors, tree, point_index, r)
    return indices(neighbors), distances(neighbors)
end

function rangesearch_sq{T}(tree::VPTree, point_index::Int, r::T)
    neighbors = NeighborList{T}(tree.num_points)
    rangesearch_sq!(neighbors, tree, point_index, r)
    return indices(neighbors), distances(neighbors)
end
