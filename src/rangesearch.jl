type NeighborList{T}
    distances::Vector{T}
    indices::Vector{Int}
    length::Int
    max_length::Int

    NeighborList(max_length::Int) = new(Vector{T}(max_length), Vector{Int}(max_length), 0, max_length)
end

clear!(neighbors::NeighborList) = neighbors.length = 0

# Is sqrt(x1) + sqrt(x2) >= sqrt(y)?
@inline function sumsqrtgt{T}(x1::T, x2::T, y::T)
    d = y - x1 - x2
    d <= 0 || d^2 <= 4*x1*x2
end

@inline function add_neighbor!{T}(neighbors::NeighborList{T}, d::T, i::Int)
    li = neighbors.length + 1
    @inbounds neighbors.distances[li] = d
    @inbounds neighbors.indices[li] = i
    neighbors.length = li
end

function add_neighbors!{T}(distance, neighbors::NeighborList{T}, point, points, node::Node{T})
    i = node.index
    i == 0 && return

    d = distance(points, i, point)::T

    add_neighbor!(neighbors, d, i)

    node.isleaf && return

    add_neighbors!(distance, neighbors, point, points, node.inside)
    add_neighbors!(distance, neighbors, point, points, node.outside)
end

function inrange!{T}(distance::Function, neighbors::NeighborList{T}, r::T, point, points, node::Node{T})
    i = node.index
    i == 0 && return

    d = distance(points, i, point)::T

    d <= r && add_neighbor!(neighbors, d, i)

    node.isleaf && return

    μ = node.radius
    # d + μ <= r
    if !sumsqrtgt(d, μ, r)
        add_neighbors!(distance, neighbors, point, points, node.inside)
    # d - μ <= r
    elseif sumsqrtgt(r, μ, d)
        inrange!(distance, neighbors, r, point, points, node.inside)
    end

    # μ - d <= r
    if sumsqrtgt(r, d, μ)
        inrange!(distance, neighbors, r, point, points, node.outside)
    end
end

function inrange!{T}(distance::Function, neighbors::NeighborList{T}, tree::VPTree, r::T, point)
    clear!(neighbors)
    inrange!(distance, neighbors, r, point, tree.points, tree.root)
end

function inrange{T}(distance::Function, tree::VPTree, r::T, point)
    neighbors = NeighborList{T}(tree.num_points)
    inrange!(distance, neighbors, r, point, tree.points, tree.root)
end
