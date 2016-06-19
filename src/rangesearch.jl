import Base: length, collect

type NeighborList{T}
    indices::Vector{Int}
    distances::Vector{T}
    length::Int
    max_length::Int

    NeighborList(max_length::Int) = new(Vector{Int}(max_length), Vector{T}(max_length), 0, max_length)
end

length(neighbors::NeighborList) = neighbors.length
clear!(neighbors::NeighborList) = neighbors.length = 0
indices(neighbors::NeighborList) = neighbors.indices[1:neighbors.length]
distances(neighbors::NeighborList) = neighbors.distances[1:neighbors.length]

# Is sqrt(x1) + sqrt(x2) >= sqrt(y)?
@inline function sumsqrtgt{T}(x1::T, x2::T, y::T)
    d = y - x1 - x2
    d <= 0 || d^2 <= 4*x1*x2
end

@inline function push_neighbor!{T}(neighbors::NeighborList{T}, index::Int, distance::T)
    i = neighbors.length + 1
    @inbounds neighbors.indices[i] = index
    @inbounds neighbors.distances[i] = distance
    neighbors.length = i
end

function push_neighbors!{T}(metric, neighbors::NeighborList{T}, node::Node{T}, point_index::Int)
    node.index == 0 && return
    node.index == point_index && return

    d = metric(point_index, node.index)::T
    push_neighbor!(neighbors, node.index, d)

    node.isleaf && return

    push_neighbors!(metric, neighbors, node.inside, point_index)
    push_neighbors!(metric, neighbors, node.outside, point_index)
end

function rangesearch!{T}(metric::Function, neighbors::NeighborList{T}, node::Node{T}, point_index::Int, r::T)
    node.index == 0 && return
    node.index == point_index && return

    d = metric(point_index, node.index)::T
    d <= r && push_neighbor!(neighbors, node.index, d)

    node.isleaf && return

    μ = node.radius
    if d + μ <= r
        push_neighbors!(metric, neighbors, node.inside, point_index)
    elseif d - μ <= r
        rangesearch!(metric, neighbors, node.inside, point_index, r)
    end

    if μ - d <= r
        rangesearch!(metric, neighbors, node.outside, point_index, r)
    end
end

function rangesearch_sq!{T}(metric::Function, neighbors::NeighborList{T}, node::Node{T}, point_index::Int, r::T)
    node.index == 0 && return
    node.index == point_index && return

    d = metric(point_index, node.index)::T
    d <= r && push_neighbor!(neighbors, node.index, d)

    node.isleaf && return

    μ = node.radius
    # d + μ <= r
    if !sumsqrtgt(d, μ, r)
        push_neighbors!(metric, neighbors, node.inside, point_index)
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
