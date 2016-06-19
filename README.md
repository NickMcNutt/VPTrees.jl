# VPTrees.jl
=====

[![Build Status](https://travis-ci.org/NickMcNutt/VPTrees.jl.svg?branch=master)](https://travis-ci.org/NickMcNutt/VPTrees.jl)
[![Coverage Status](https://coveralls.io/repos/NickMcNutt/VPTrees.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/NickMcNutt/VPTrees.jl?branch=master)
[![codecov.io](http://codecov.io/github/NickMcNutt/VPTrees.jl/coverage.svg?branch=master)](http://codecov.io/github/NickMcNutt/VPTrees.jl?branch=master)

## About

### Overview

* A Julia package for creating and searching Vantage Point Trees

* Functionality includes k-nearest neighbors (kNN) search and range search

* Code is focused on high performance and high generality

* Licensed under the [MIT License](https://opensource.org/licenses/MIT)

### What are VP Trees?

Vantage Point Trees (VP Trees) are similar to KD trees; however, constructing a KD tree requires knowledge of the structure of the space in which the data points live.  VP Trees on the other hand require nothing more than a [metric](https://en.wikipedia.org/wiki/Metric_space) that gives the distance between two points.

### How to install

In Julia:
```julia
Pkg.clone("https://github.com/NickMcNutt/VPTrees.jl")
```

## Examples

```julia

# Create some points in 3D space
n = 100000
points = randn(3, n)
point_indices = collect(1:n)

# Define a metric (Euclidean distance in this case)
function metric(i::Int, j::Int)
    x = points[1, i] - points[1, j]
    y = points[2, i] - points[2, j]
    z = points[3, i] - points[3, j]

    return sqrt(x^2 + y^2 + z^2)
end

# Construct a VP Tree containing all of the points
tree = VPTree(metric, point_indices)

# Find the 10 nearest points to point 123
indices, distances = knn(tree, 123, 10)

# Find all points within a distance of 1.5 from point 123
indices, distances = rangesearch(tree, 123, 1.5)
```

## Functionality

### Typical usage

```julia
VPTree(metric::Function, point_indices::Vector{Int})
```

Construct a VP Tree from a set of points indexed by `point_indices`.  The user-supplied function `metric(i::Int, j::Int)` takes the indices of two points and returns the distance between them.

```julia
knn(tree::VPTree, point_index::Int, k::Int)
```

Search a VP Tree for the `k` nearest neighbors to the point indexed by `point_index`.  Returns a list of neighbor indices and a list of their corresponding distances.

```julia
rangesearch(tree:VPTree, point_index::Int, radius::T)
```

Search a VP Tree for all points within a distance of `radius` from the point indexed by `point_index`.  Returns a list of neighbor indices and a list of their corresponding distances.


### High performance usage

#### Preallocated memory

If a VP Tree is intended to be used for a significant number of search queries, then one neighbor list can be used for all searches in order to reduce the amount of time spent allocating/deallocating memory.

```julia
knn!(neighbors::NeighborList, tree::VPTree, point_index::Int, num_neighbors::Int)

rangesearch!(neighbors::NeighborList, tree:VPTree, point_index::Int, radius::T)
```

These functions are the same as `knn()` and `rangesearch()`, but use the neighbor list `neighbors` to store the results of the search query.  The initial contents of the list are cleared.

```julia
NeighborList{T}(num_neighbors::Int)
```

Construct an empty neighbor list with room for up to `num_neighbors` neighbors.  The neighbor distances are of type `T`.

```julia
indices(neighbors::NeighborList)
```

Return a list of neighbor indices.

```julia
distances(neighbors::NeighborList)
```

Return a list of neighbor distances.

#### Squared distances

```julia
knn_sq(), knn_sq!(), rangesearch_sq(), rangesearch_sq!()
```

In some cases, the performance bottleneck for creating and searching VP Trees is the `sqrt()` function.  It is often faster to work directly with the squared distances between points instead.  These variants of the search functions assume that all distance values are actually *squared* distance values.

### Additional functionality
```julia
knn!(metric::Function, neighbors::NeighborList, node::Node, point_index::Int, num_neighbors::Int)

rangesearch!(metric::Function, neighbors::NeighborList, node::Node, point_index::Int, r::T)
```
Searches a given `node` using the metric function `metric` and neighbor list `neighbors`.


## Current status

### Complete

* VP Tree construction
* k-nearest neighbor search for points that exist within the VP Tree
* Range search for points that exist within the VP Tree

### Todo

* k-nearest neighbor search for points not in the VP Tree
* Range search for points not in the VP Tree
* Insert, replace, and delete points from an existing VP Tree
* Figure out why `knn!()` and `rangesearch!()` allocate memory.  These functions should not be allocating any memory.
* Add support for using the `Distances.jl` package in order to streamline the usage of common metrics.
