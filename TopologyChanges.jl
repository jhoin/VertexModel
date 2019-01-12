module TopologyChanges

include("dcel.jl")
using .DCEL

# Perform t1 transition
# Arguments:edge object
# Return: edge object updated
function t1transition!(focal_edge)

    # Perform the topological changes for the focal edge
    t1topology!(focal_edge)

    # Perform the topological changes for the focal twin
    t1topology!(focal_edge.twinEdge)

    # Rotate points
    p1 = focal_edge.originVertex
    p2 = focal_edge.nextEdge.originVertex
end # t1transition

# Make the topological changes for the t1 transition
# Arguments:edge object
# Return: edge object updated
function t1topology!(focal_edge)

    # Check if edge is not incident edge of cell
    focal_cell = focal_edge.containCell
    if focal_cell.incEdge == focal_edge
        focal_edge.incEdge = focal_edge.nextEdge
    end

    # Edge operations
    edge = focal_edge.prevEdge
    edgeremove!(focal_edge)
    addedgeat!(focal_edge, edge)
end # t1topolgy

# Remove edge from cell
# Arguments: edge object
# Return: edge object
function edgeremove!(edge)
    edge.prevEdge.nextEdge = edge.nextEdge
end # edgeremove

# Add edge into cell before indicated edge
# Arguments: edge to be added, and edge that will come after the given edge
# Return: edge object
function addedgeat!(focal_edge,next_edge)
    next_edge.originVertex = focal_edge.nextEdge.originVertex
    focal_edge.nextEdge = next_edge
    focal_edge.prevEdge = next_edge.prevEdge
    next_edge.prevEdge = focal_edge
end # addedgeat

end # module
