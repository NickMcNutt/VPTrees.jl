module VantagePointTrees

using Distances

export
    # construct.jl
    PeriodicSqEuclidean,
    evaluate,

    Node,
    VPTree,
    partition!,
    quickselect!,
    distances!,
    #workmem

    # query.jl
    knn!,

    # display.jl
    printnode,
    show,
    display,
    writemime

include("construct.jl")
include("query.jl")
include("display.jl")

end
