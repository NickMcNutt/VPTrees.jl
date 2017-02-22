import Base: show, display, length, collect, indices

type NeighborList{T}
    indices::Vector{Int}
    distances::Vector{T}
    length::Int
    max_length::Int

    NeighborList{T}(max_length::Int) where T = new(Vector{Int}(max_length), Vector{T}(max_length), 0, max_length)
end

length(neighbors::NeighborList) = neighbors.length
indices(neighbors::NeighborList) = neighbors.indices[1:neighbors.length]
distances(neighbors::NeighborList) = neighbors.distances[1:neighbors.length]

function show(io::IO, neighbors::NeighborList)
    for i in 1:neighbors.length
        @printf(io, "%-10d%.4f\n", neighbors.indices[i], neighbors.distances[i])
    end

    for i in 1:neighbors.max_length-neighbors.length
        @printf(io, "*\n")
    end
end

display(neighbors::NeighborList) = show(neighbors)

clear!(neighbors::NeighborList) = neighbors.length = 0

function insert_neighbor!{T}(neighbors::NeighborList{T}, index::Int, distance::T) 
    ni = neighbors.indices
    nd = neighbors.distances
    k = neighbors.length

    if k < neighbors.max_length
        neighbors.length += 1
        k += 1
        @inbounds nd[k] = T(Inf)
    end

    j = 1
    @inbounds while distance > nd[j]
        j += 1
    end

    for m in k:-1:j+1
        @inbounds nd[m] = nd[m - 1]
        @inbounds ni[m] = ni[m - 1]
    end

    @inbounds ni[j] = index
    @inbounds nd[j] = distance

    @inbounds τ = k == neighbors.max_length ? nd[k] : T(Inf)

    return τ
end

@inline function add_neighbor!{T}(neighbors::NeighborList{T}, index::Int, distance::T)
    k = neighbors.length + 1
    @inbounds neighbors.indices[k] = index
    @inbounds neighbors.distances[k] = distance
    neighbors.length = k
end

function add_neighbors!{T}(metric, neighbors::NeighborList{T}, node::Node{T}, point_index::Int)
    node.index == 0 && return

    #if node.index != point_index
    #if d > 0.001

    d = metric(node.index, point_index)::T
    add_neighbor!(neighbors, node.index, d)

    node.isleaf && return

    add_neighbors!(metric, neighbors, node.inside, point_index)
    add_neighbors!(metric, neighbors, node.outside, point_index)
end

# Is sqrt(x1) + sqrt(x2) >= sqrt(y)?
@inline function sumsqrtgt{T}(x1::T, x2::T, y::T)
    d = y - x1 - x2
    d <= 0 || d^2 <= 4*x1*x2
end
