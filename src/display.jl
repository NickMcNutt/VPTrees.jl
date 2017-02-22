import Base.show, Base.display

function printnode(io::IO, node::Node, maxdepth::Int = 4, depth::Int = 0)
    junction = ifelse(depth > 0, "├", "┌")
    tab = "   │   " ^ depth * "   $(junction)── "
    if node.isleaf
        if node.index != 0
            @printf(io, "%s(%d)\n", tab, node.index)
        else
            @printf(io, "%s*\n", tab)
        end
    elseif depth == maxdepth
        println(io, tab, "...")
    else
        @printf(io, "%s%d, %#-.1f\n", tab, node.index, node.radius)
        #@printf(io, "%s[%d, %d], %d, %#-.1f\n", tab, node.i, node.j, node.index, sqrt(node.radius))
        printnode(io, node.inside, maxdepth, depth + 1)
        printnode(io, node.outside, maxdepth, depth + 1)
    end
end

show(io::IO, node::Node, args...) = printnode(io, node, args...)
display(node::Node) = show(node)

show(io::IO, vptree::VPTree, args...) = show(io, vptree.root, args...)
display(vptree::VPTree) = show(vptree.root)
