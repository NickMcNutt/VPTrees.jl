using VPTrees, DataStructures
using Base.Test

function euclidean(coords::Matrix, i::Int, j::Int)
    @inbounds x = coords[1, i] - coords[1, j]
    @inbounds y = coords[2, i] - coords[2, j]
    @inbounds z = coords[3, i] - coords[3, j]

    return sqrt(x^2 + y^2 + z^2)
end

function euclidean_periodic(coords::Matrix, i::Int, j::Int)
    xbw, ybw, zbw = 1.0, 1.0, 1.0

    xhbw = xbw / 2
    yhbw = ybw / 2
    zhbw = zbw / 2

    @inbounds x = coords[1, i] - coords[1, j]
    @inbounds y = coords[2, i] - coords[2, j]
    @inbounds z = coords[3, i] - coords[3, j]

    if x < -xhbw x += xbw elseif x > xhbw x -= xbw end
    if y < -yhbw y += ybw elseif y > yhbw y -= ybw end
    if z < -zhbw z += zbw elseif z > zhbw z -= zbw end

    return sqrt(x^2 + y^2 + z^2)
end

function is_consistent{T}(metric::Function, node::Node{T}, point_index::Int, τ::T, isless)
    node.index == 0 && return true

    μ = metric(node.index, point_index)
    #println("$(point_index) and $(node.index): $μ <= $τ is $(μ <= τ) and isless = $isless")
    ((μ <= τ) ⊻ isless) && return false

    node.isleaf && return true

    is_consistent(metric, node.inside, point_index, τ, isless) || return false
    is_consistent(metric, node.outside, point_index, τ, isless) || return false

    return true
end

function is_consistent{T}(metric::Function, node::Node{T})
    node.isleaf && return true

    is_consistent(metric, node.inside, node.index, node.radius, true) || return false
    is_consistent(metric, node.outside, node.index, node.radius, false) || return false
    #println()
    is_consistent(metric, node.inside) || return false
    is_consistent(metric, node.outside) || return false

    return true
end

is_consistent(tree::VPTree) = is_consistent(tree.metric, tree.root)

create_metric(distance::Function, points) = (i::Int, j::Int) -> distance(points, i, j)

function test_knn(distance::Function, num_points::Int, num_iters::Int)
    for s in 1:num_iters
        srand(s)
        points = rand(3, num_points)
        point_indices = collect(1:num_points)
        metric = create_metric(distance, points)

        num_neighbors = 5

        tree = VPTree(metric, point_indices)
        neighbors = NeighborList{Float64}(num_neighbors)

        for i in 1:num_points
            knn!(neighbors, tree, i)
            indices_test = indices(neighbors)

            pq = PriorityQueue(Int, Float64, Base.Order.Reverse)
            for j in 1:num_points
                i == j && continue

                r = metric(i, j)

                if length(pq) < num_neighbors
                    enqueue!(pq, j, r)
                elseif r < last(peek(pq))
                    dequeue!(pq)
                    enqueue!(pq, j, r)
                end
            end

            indices_actual = keys(pq)

            isempty(symdiff(indices_test, indices_actual)) || return false
        end
    end

    return true
end

function test_rangesearch(distance::Function, num_points::Int, num_iters::Int)
    for s in 1:num_iters
        srand(s)
        points = rand(3, num_points)
        point_indices = collect(1:num_points)
        metric = create_metric(distance, points)

        tree = VPTree(metric, point_indices)
        neighbors = NeighborList{Float64}(num_points)

        r = rand()
        for i in 1:num_points
            rangesearch!(neighbors, tree, i, r)
            indices_test = indices(neighbors)

            indices_actual = filter(j -> metric(i, j) <= r, 1:num_points)

            isempty(symdiff(indices_test, indices_actual)) || return false
        end
    end

    return true
end

function test_integrity(distance::Function, num_points::Int, num_iters::Int)
    for s in 1:num_iters
        srand(s)
        points = rand(3, num_points)
        point_indices = collect(1:num_points)
        metric = create_metric(distance, points)

        tree = VPTree(metric, point_indices)
        is_consistent(tree) || return false
    end

    return true
end

@testset "VP Tree tests" begin
    const num_points = 1000
    const num_iters = 5

    @testset "VP Tree integrity test" begin
        @test test_integrity(euclidean, num_points, num_iters)
        @test test_integrity(euclidean_periodic, num_points, num_iters)
    end

    @testset "k-nearest neighbors test" begin
        @test test_knn(euclidean, num_points, num_iters)
        @test test_knn(euclidean_periodic, num_points, num_iters)
    end

    @testset "Range search test" begin
        @test test_rangesearch(euclidean, num_points, num_iters)
        @test test_rangesearch(euclidean_periodic, num_points, num_iters)
    end
end
