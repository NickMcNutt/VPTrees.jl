module VPTrees

export
    # construct.jl
    Node,
    VPTree,

    # neighborlist.jl
    NeighborList,
    length,
    indices,
    distances,

    # knn.jl
    knn!,
    knn_sq!,

    # rangesearch.jl
    rangesearch!,
    rangesearch,
    rangesearch_sq!,
    rangesearch_sq,

    # display.jl
    printnode,
    show,
    display

include("construct.jl")
include("neighborlist.jl")
include("knn.jl")
include("rangesearch.jl")
include("display.jl")

end
