module VPTrees

export
    # construct.jl
    Node,
    VPTree,

    # knn.jl
    knn!,

    # rangesearch.jl
    NeighborList,
    inrange!,
    inrange,

    # display.jl
    printnode,
    show,
    display

include("construct.jl")
include("knn.jl")
include("rangesearch.jl")
include("display.jl")

end
