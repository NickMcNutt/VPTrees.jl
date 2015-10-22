import Base.show, Base.display, Base.writemime

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
        @printf(io, "%s, %d, %#-.1f\n", tab, node.index, sqrt(node.distance))
        #@printf(io, "%s[%d, %d], %d, %#-.1f\n", tab, node.i, node.j, node.index, sqrt(node.distance))
        printnode(io, node.left, maxdepth, depth + 1)
        printnode(io, node.right, maxdepth, depth + 1)
    end
end

show(io::IO, node::Node, args...) = printnode(io, node, args...)
display(node::Node) = show(node)
writemime(io::IO, ::MIME"text/plain", node::Node, args...) = printnode(io, node, args...)

show(io::IO, vptree::VPTree, args...) = show(io, vptree.root, args...)
display(vptree::VPTree) = show(vptree.root)
writemime(io::IO, ::MIME"text/plain", vptree::VPTree, args...) = writemime(io, MIME"text/plain"(), vptree.root, args...)
