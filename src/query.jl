# Is sqrt(x1) + sqrt(x2) >= sqrt(y)?
function sumsqrtgt{T<:Real}(x1::T, x2::T, y::T)
    d = y - x1 - x2
    if d > 0 && d^2 > 4*x1*x2
        return false
    else
        return true
    end
end

function initneighbors!{T<:Real}(nd::AbstractVector{T}, k::Int)
    @inbounds @simd for i in 1:k
        nd[i] = Inf
    end
end

function updateneighbors!{T<:Real}(nd::AbstractVector{T}, ni::AbstractVector{Int}, d::T, i::Int, k::Int)
    j = 1
    @inbounds while d > nd[j]
        j += 1
    end

    @inbounds for m in k:-1:j+1
        nd[m] = nd[m - 1]
        ni[m] = ni[m - 1]
    end

    @inbounds nd[j] = d
    @inbounds ni[j] = i

    @inbounds τ = nd[k]
    return τ
end

function knn!{T<:Real}(nd::AbstractVector{T}, ni::AbstractVector{Int}, τ::T, metric::Metric, points, node::Node, k::Int, x::T, y::T, z::T)
    @inbounds begin
        i = node.index
        i == 0 && return τ

        d::T = evaluate(metric, points, i, x, y, z)::T

        if d < τ
            τ = updateneighbors!(nd, ni, d, i, k)
        end

        node.isleaf && return τ

        μ = node.distance
        # if τ + d <= μ
        if !sumsqrtgt(τ, d, μ)
            τ = knn!(nd, ni, τ, metric, points, node.left, k, x, y, z)
        else
            τ = knn!(nd, ni, τ, metric, points, node.right, k, x, y, z)
            # if d <= τ + μ
            if sumsqrtgt(τ, μ, d)
                τ = knn!(nd, ni, τ, metric, points, node.left, k, x, y, z)
            end
        end

        return τ
    end
end

function knn!{T<:Real}(nd::AbstractVector{T}, ni::AbstractVector{Int}, tree::VPTree, k::Int, x::T, y::T, z::T)
    initneighbors!(nd, k)
    knn!(nd, ni, T(Inf), tree.metric, tree.points, tree.root, k, x, y, z)
end
