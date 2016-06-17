# Is sqrt(x1) + sqrt(x2) >= sqrt(y)?
@inline function sumsqrtgt{T}(x1::T, x2::T, y::T)
    d = y - x1 - x2
    d <= 0 || d^2 <= 4*x1*x2
end

function initneighbors!{T}(nd::AbstractVector{T}, k::Int)
    @inbounds @simd for i in 1:k
        nd[i] = Inf
    end
end

function updateneighbors!{T}(nd::AbstractVector{T}, ni::AbstractVector{Int}, d::T, i::Int, k::Int)
    j = 1
    @inbounds while d > nd[j]
        j += 1
    end

    for m in k:-1:j+1
        @inbounds nd[m] = nd[m - 1]
        @inbounds ni[m] = ni[m - 1]
    end

    @inbounds nd[j] = d
    @inbounds ni[j] = i

    @inbounds τ = nd[k]
    return τ
end

function knn!{T}(distance::Function, nd::AbstractVector{T}, ni::AbstractVector{Int}, τ::T, points, node::Node{T}, k::Int)
    i = node.index
    i == 0 && return τ

    d = distance(points, i)::T

    if d < τ
        τ = updateneighbors!(nd, ni, d, i, k)
    end

    node.isleaf && return τ

    μ = node.radius
    # τ + d <= μ
    if !sumsqrtgt(τ, d, μ)
        τ = knn!(distance, nd, ni, τ, points, node.inside, k)
    else
        τ = knn!(distance, nd, ni, τ, points, node.outside, k)
        # d <= τ + μ
        if sumsqrtgt(τ, μ, d)
            τ = knn!(distance, nd, ni, τ, points, node.inside, k)
        end
    end

    return τ
end

function knn!{T}(distance::Function, nd::AbstractVector{T}, ni::AbstractVector{Int}, tree::VPTree, k::Int)
    initneighbors!(nd, k)
    knn!(distance, nd, ni, T(Inf), tree.points, tree.root, k)
end

function addneighbors!{T}(distance, nd::AbstractVector{T}, ni::AbstractVector{Int}, points, node::Node{T})
    i = node.index
    i == 0 && return

    d = distance(points, i)::T

    push!(nd, d)
    push!(ni, i)

    node.isleaf && return

    addneighbors!(distance, nd, ni, points, node.inside)
    addneighbors!(distance, nd, ni, points, node.outside)
end

function inrange!{T}(distance::Function, nd::AbstractVector{T}, ni::AbstractVector{Int}, r::T, points, node::Node{T})
    i = node.index
    i == 0 && return

    d = distance(points, i)::T

    if d <= r
        push!(nd, d)
        push!(ni, i)
    end

    node.isleaf && return

    μ = node.radius
    # d + μ <= r
    if !sumsqrtgt(d, μ, r)
        addneighbors!(distance, nd, ni, points, node.inside)
    # d - μ <= r
    elseif sumsqrtgt(r, μ, d)
        inrange!(distance, nd, ni, r, points, node.inside)
    end

    # μ - d <= r
    if sumsqrtgt(r, d, μ)
        inrange!(distance, nd, ni, r, points, node.outside)
    end
end

function inrange{T}(distance::Function, tree::VPTree, r::T)
    nd = Vector{T}()
    ni = Vector{Int}()
    inrange!(distance, nd, ni, r, tree.points, tree.root)
    return nd, ni
end
