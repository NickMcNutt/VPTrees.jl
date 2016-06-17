module VPTrees

export
    # construct.jl
    Node,
    VPTree,

    # query.jl
    knn!,
    inrange!,
    inrange,

    # display.jl
    printnode,
    show,
    display

include("construct.jl")
include("query.jl")
include("display.jl")

end
