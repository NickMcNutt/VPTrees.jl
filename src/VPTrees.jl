module VPTrees

export
    # construct.jl
    Node,
    VPTree,

    # knn.jl
    knn!,
    knn_sq!,

    # rangesearch.jl
    NeighborList,
    rangesearch!,
    rangesearch,
    rangesearch_sq!,
    rangesearch_sq,

    # display.jl
    printnode,
    show,
    display

include("construct.jl")
include("knn.jl")
include("rangesearch.jl")
include("display.jl")

end
