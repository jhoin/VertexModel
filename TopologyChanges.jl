module TopologyChanges

using ..DCEL
using ..PointGeometry

# Check for topological changes in mesh
# Arguments: mesh
# Return: updated mesh
function update_topology!(mesh::Dcel, minlen::Float64)

    # Check for t1 transitions
    for i in 1:size(mesh.listEdge,1)
        if mesh.listEdge[i].border
            continue
        end
        if mesh.listEdge[i].edgeLen < minlen
            t1transition!(mesh.listEdge[i], minlen)
            #cell = mesh.listEdge[i].containCell
            #cell_division!(mesh, cell)
            #break
        end
    end

    #for i in 1:size(mesh.listCell,1)
    #    cell_division!(mesh, mesh.listCell[i])
    #end


end #update_topology
export update_topology!

# Perform t1 transition
# Arguments:edge object
# Return: edge object updated
function t1transition!(focal_edge::Hedge, minlen::Float64)

    # Assert if neighbor edges are in the border

    # Rotate points
    p1 = focal_edge.originVertex
    p2 = focal_edge.nextEdge.originVertex
    rotatepoints!(p1, p2)

    # Perform the topological changes for the focal edge
    t1topology!(focal_edge)

    # Perform the topological changes for the focal twin
    t1topology!(focal_edge.twinEdge)
    if focal_edge.prevEdge.border println("here") end

    # Update edge length
    extend_edge(focal_edge, minlen)
end # t1transition
export t1transition!


# Make the topological changes for the t1 transition
# Arguments:edge object
# Return: edge object updated
function t1topology!(focal_edge::Hedge)

    # Check if edge is not incident edge of cell
    focal_cell = focal_edge.containCell
    if focal_cell.incEdge == focal_edge
        focal_cell.incEdge = focal_edge.nextEdge
    end
    println(isdefined(focal_edge, :prevEdge))

    # Edge operations
    edge = focal_edge.prevEdge.twinEdge
    edgeremove!(focal_edge, focal_edge.prevEdge)
    addedgeat!(focal_edge, edge)
end # t1topolgy

# Remove edge from cell
# Arguments: edge object
# Return: edge object
function edgeremove!(edge::Hedge, prevtwin::Hedge)
    edge.prevEdge.nextEdge = edge.nextEdge
    edge.containCell = prevtwin.twinEdge.containCell
end # edgeremove

# Add edge into cell before indicated edge
# Arguments: edge to be added, and edge that will come after the given edge
# Return: edge object
function addedgeat!(focal_edge::Hedge, next_edge::Hedge)
    next_edge.originVertex = focal_edge.nextEdge.originVertex
    focal_edge.nextEdge = next_edge
    focal_edge.prevEdge = next_edge.prevEdge
    next_edge.prevEdge.nextEdge = focal_edge
    next_edge.prevEdge = focal_edge
end # addedgeat

# create an edge cell after the edge
# Arguments: mesh, and edge that will come before the given edge and origin
# Return: mesh object with new
function addedgeat!(mesh::Dcel, prev_edge::Hedge, vert::Vertex)
    edge = Hedge() # Edge created by the intersection
    push!(mesh.listEdge, edge)
    edge.originVertex = vert
    edge.border = prev_edge.border
    edge.prevEdge = prev_edge
    edge.nextEdge = prev_edge.nextEdge

    prev_edge.nextEdge.prevEdge = edge
    prev_edge.nextEdge = edge
    edge.containCell = prev_edge.containCell
    return edge
end # addedgeat!

# create an edge cell between two edges
# Arguments: mesh, previous and next edge and the origin vertex
# Return: mesh object with new
function addedgeat!(mesh::Dcel, prev_edge::Hedge, next_edge::Hedge, vert::Vertex)
    edge = Hedge()
    push!(mesh.listEdge, edge)
    edge.originVertex = vert
    edge.prevEdge = prev_edge
    prev_edge.nextEdge = edge
    edge.nextEdge = next_edge
    next_edge.prevEdge = edge
    edge.containCell = prev_edge.containCell
    return edge
end # addedgeat!

# Perform a cell division
# Arguments: mesh object
# Return: mesh object updated
function cell_division!(mesh::Dcel, cell::Cell)

    # Get the shortest axis, aka, the vertex in the shortest axis
    c2 = shortest_axis(cell)

    # loop over edges of cell and look for the points where the axis crosses then edge
    edges_crossed = Vector{Hedge}()
    verts_crossing = Vector{Vertex}()
    for edge in cell
        inter = intersection(edge.originVertex, edge.nextEdge.originVertex, cell.centroid, c2)
        if !isnan(inter.x) && is_between(edge.originVertex, edge.nextEdge.originVertex, inter)
            push!(verts_crossing, inter)
            push!(mesh.listVert, inter)
            push!(edges_crossed, edge)
        end
    end

    # Topological changes in edge
    split1 = addedgeat!(mesh, edges_crossed[1], verts_crossing[1])
    split2 = addedgeat!(mesh, edges_crossed[2], verts_crossing[2])
    if !edges_crossed[1].border
        split1_twin = addedgeat!(mesh, edges_crossed[1].twinEdge.prevEdge, edges_crossed[1].twinEdge.originVertex)
        edges_crossed[1].twinEdge.originVertex = verts_crossing[1]
        split1_twin.twinEdge = split1
        split1.twinEdge = split1_twin
        #split1_twin.twinEdge, edges_crossed[1].twinEdge = edges_crossed[1], split1_twin
    end

    if !edges_crossed[2].border
        split2_twin = addedgeat!(mesh, edges_crossed[2].twinEdge.prevEdge, edges_crossed[2].twinEdge.originVertex)
        edges_crossed[2].twinEdge.originVertex = verts_crossing[2]
        split2_twin.twinEdge = split2
        split2.twinEdge = split2_twin
    end

    focal1 = addedgeat!(mesh, edges_crossed[1], split2, verts_crossing[1])
    focal2 = addedgeat!(mesh, edges_crossed[2], split1, verts_crossing[2])
    focal1.twinEdge, focal2.twinEdge = focal2, focal1

    for i in 1:size(mesh.listEdge, 1)
        newedgelen!(mesh.listEdge[i])
    end

    # Topological changes in cells
    new_cell = Cell()
    push!(mesh.listCell, new_cell)
    new_cell.incEdge = split1
    focal2.containCell = new_cell
    for edge in new_cell
       edge.containCell = new_cell
    end

end # cell_division!
export cell_division!


end # module
