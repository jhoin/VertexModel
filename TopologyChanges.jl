module TopologyChanges

using ..Meshes
using ..PointGeometry

# Check for topological changes in mesh
# Arguments: mesh
# Return: updated mesh
function update_topology!(mesh::Mesh, minlen::Float64, K::Float64,minarea::Float64)

    for i in 1:size(mesh.cells,1)
        cell = mesh.cells[i]
        if cell.areaCell > 2.5*K
            cell_division!(mesh, cell)
        end
        if cell.areaCell < minarea
            halt = check_t2(cell)
            halt && continue
            length(cell) < 4 && cell_extrusion!(mesh, cell)
        end
    end

    # Check for t1 transitions
    for i in 1:size(mesh.edges,1)
        halt = check_t1(mesh.edges[i])
        halt && continue
        if mesh.edges[i].edgeLen < minlen
            t1transition!(mesh.edges[i], minlen)
            break
        end
    end

end #update_topology
export update_topology!

# Perform t1 transition
# Arguments:edge object
# Return: edge object updated
function t1transition!(focal_edge::Hedge, minlen::Float64)

    # Rotate points
    p1 = focal_edge.originVertex
    p2 = focal_edge.nextEdge.originVertex
    rotatepoints!(p1, p2)

    # Store the affected cells here
    aff_cells = Vector{Cell}(undef,4)

    # Perform the topological changes for the focal edge
    aff_cells[1:2] = t1topology!(focal_edge)

    # Perform the topological changes for the focal twin
    aff_cells[3:4] = t1topology!(focal_edge.twinEdge)

    # Update the edges leaving the vertices involved
    setleavingedge!(focal_edge.nextEdge.originVertex, focal_edge.twinEdge.nextEdge, focal_edge.nextEdge)
    focal_twin= focal_edge.twinEdge
    setleavingedge!(focal_twin.nextEdge.originVertex, focal_twin.twinEdge.nextEdge, focal_twin.nextEdge)

    # Update edge length
    extend_edge(focal_edge, minlen)

    # Update cell measures
    for i in 1:length(aff_cells)
        updatecell!(aff_cells[i])
    end
end # t1transition
export t1transition!

function check_t1(edge::Hedge)
    cannot_perform = false
    if edge.border cannot_perform = true end
    if edge.nextEdge.border cannot_perform = true end
    if edge.prevEdge.border  cannot_perform = true end
    return cannot_perform
end

function check_t2(cell::Cell)
    halt = false
    for edge in cell
        if edge.border halt = true end
    end
    return halt
end
# Make the topological changes for the t1 transition
# Arguments:edge object
# Return: edge object updated
function t1topology!(focal_edge::Hedge)

    # Check if edge is not incident edge of cell
    focal_cell = focal_edge.containCell
    if focal_cell.incEdge == focal_edge
        focal_cell.incEdge = focal_edge.nextEdge
    end

    # Edge operations
    edge = focal_edge.prevEdge.twinEdge
    adj_cell = edge.containCell
    #setleavingedge!(focal_edge.nextEdge.originVertex, focal_edge.twinEdge.nextEdge, focal_edge.nextEdge)
    edgeremove!(focal_edge, edge)
    addedgeat!(focal_edge, edge)
    return [focal_cell, adj_cell]
end # t1topolgy

# Remove edge from cell
# Arguments: edge object
# Return: edge object
function edgeremove!(edge::Hedge, prevtwin::Hedge)
    edge.containCell = prevtwin.containCell
    edge.prevEdge.nextEdge = edge.nextEdge
    edge.nextEdge.prevEdge = edge.prevEdge
end # edgeremove

function edgeremove!(edge::Hedge)
    edge.prevEdge.nextEdge = edge.nextEdge
    edge.nextEdge.prevEdge = edge.prevEdge
end

# Add edge into cell before indicated edge
# Arguments: edge to be added, and edge that will come after the given edge
# Return: edge object
function addedgeat!(focal_edge::Hedge, next_edge::Hedge)
    next_edge.originVertex = focal_edge.nextEdge.originVertex
    #setleavingedge!(next_edge.originVertex, focal_edge.twinEdge.prevEdge, next_edge)
    focal_edge.nextEdge = next_edge
    focal_edge.prevEdge = next_edge.prevEdge
    next_edge.prevEdge.nextEdge = focal_edge
    next_edge.prevEdge = focal_edge
end # addedgeat

# create an edge cell after the edge
# Arguments: mesh, and edge that will come before the given edge and origin
# Return: mesh object with new
function addedgeat!(mesh::Mesh, prev_edge::Hedge, vert::Vertex)
    edge = Hedge() # Edge created by the intersection
    #addtomesh!(mesh, edge)
    edge.originVertex = vert
    setleavingedge!(vert, edge)
    edge.border = prev_edge.border
    edge.prevEdge = prev_edge
    edge.nextEdge = prev_edge.nextEdge

    prev_edge.nextEdge.prevEdge = edge
    prev_edge.nextEdge = edge
    edge.containCell = prev_edge.containCell
    return edge
end # addedgeat!

function addedgeat!(mesh::Mesh, prev_edge::Hedge)
    edge = Hedge() # Edge created by the intersection
    addtomesh!(mesh, edge)
    #edge.originVertex = vert
    #setleavingedge!(vert, edge)
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
function addedgeat!(mesh::Mesh, prev_edge::Hedge, next_edge::Hedge, vert::Vertex)
    edge = Hedge()
    addtomesh!(mesh, edge)
    edge.originVertex = vert
    setleavingedge!(vert, edge)
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
function cell_division!(mesh::Mesh, cell::Cell)

    # Get the shortest axis, aka, the vertex in the shortest axis
    c2 = shortest_axis(cell)

    # loop over edges of cell and look for the points where the axis crosses then edge
    edges_crossed, verts_crossing = intersectcell!(mesh, cell, c2)

    # Topological changes in edge
    affected_cells = Vector{Cell}()
    push!(affected_cells, cell)
    split1 = addedgeat!(mesh, edges_crossed[1], verts_crossing[1])
    split2 = addedgeat!(mesh, edges_crossed[2], verts_crossing[2])
    if !edges_crossed[1].border
        twin_split1= split_twin!(mesh, edges_crossed[1], split1, verts_crossing[1])
        push!(affected_cells, edges_crossed[1].twinEdge.containCell)
    end

    if !edges_crossed[2].border
        twin_split1= split_twin!(mesh, edges_crossed[2], split2, verts_crossing[2])
        push!(affected_cells, edges_crossed[2].twinEdge.containCell)
    end

    focal1 = addedgeat!(mesh, edges_crossed[1], split2, verts_crossing[1])
    focal2 = addedgeat!(mesh, edges_crossed[2], split1, verts_crossing[2])
    focal1.twinEdge, focal2.twinEdge = focal2, focal1

    for i in 1:size(mesh.edges, 1)
        newedgelen!(mesh.edges[i])
    end

    # Topological changes in cells
    new_cell = Cell()
    new_cell.centroid = Vertex([0.0,0.0])
    addtomesh!(mesh, new_cell)
    new_cell.area_elast = cell.area_elast
    new_cell.perim_elast = cell.perim_elast
    new_cell.incEdge = focal2
    focal2.containCell = new_cell
    for edge in new_cell
       edge.containCell = new_cell
    end
    push!(affected_cells, new_cell)
    for i in 1:length(affected_cells)
        updatecell!(affected_cells[i])
    end
    setleavingedges!(affected_cells)

end # cell_division!
export cell_division!

# Get the edges and verts intersecting a cell for division
function intersectcell!(mesh::Mesh, cell::Cell, vert::Vertex)
    edges_crossed = Vector{Hedge}()
    verts_crossing = Vector{Vertex}()
    no_cross = true
    for edge in cell
        inter = intersection(edge.originVertex, edge.nextEdge.originVertex, cell.centroid, vert)
        if !isnan(inter.x) && is_between(edge.originVertex, edge.nextEdge.originVertex, inter)
            no_cross = false
            push!(verts_crossing, inter)
            push!(mesh.vertices, inter)
            push!(edges_crossed, edge)
            inter.treatBoundary = Mobile(inter)
        end
    end
    return edges_crossed, verts_crossing
end

# Split the twin edge of an edge already split (cell division step)
# Arguments: mesh, edge being crossed, edge created by the crossing and the intersection vert
# Return: edge that is the twin of the
function split_twin!(mesh::Mesh, crossed::Hedge, split::Hedge, vert::Vertex)
    #split_twin = addedgeat!(mesh, crossed.twinEdge.prevEdge, vert)
    split_twin = addedgeat!(mesh, crossed.twinEdge.prevEdge)
    split_twin.originVertex = crossed.twinEdge.originVertex
    crossed.twinEdge.originVertex = vert
    setleavingedge!(crossed.twinEdge.originVertex, crossed.twinEdge, split_twin)
    setleavingedge!(vert, crossed.twinEdge)
    split_twin.twinEdge, split.twinEdge = split, split_twin
    split.twinEdge = split_twin
    return split_twin
end

function cell_extrusion!(mesh::Mesh, cell::Cell)

    for edge in cell
        edge.border && return
    end
    # pick the first vertex and move to the center
    vert = cell.incEdge.originVertex
    vert.x, vert.y = cell.centroid.x, cell.centroid.y

    # set leaving edges
    extrude_leavingedges!(cell)

    # remove twins from neighbor cells
    extrude_twins!(cell)

    delete_elements!(mesh, cell)
end

function extrude_leavingedges!(cell::Cell)
    vert = cell.incEdge.originVertex
    for i in 1:3
        edge = vert.leavingEdges[i]
        if edge.containCell == cell
            replace = edge.twinEdge.prevEdge.twinEdge
            replace.originVertex = vert
            setleavingedge!(vert, edge, replace)
        elseif edge.twinEdge.containCell == cell
            replace = edge.nextEdge
            replace.originVertex = vert
            setleavingedge!(vert, edge, replace)
        end
    end
end

function extrude_twins!(cell::Cell)
    for edge in cell
        edgeremove!(edge.twinEdge)
    end
end

function delete_elements!(mesh::Mesh, cell::Cell)
    keepvert = true
    for edge in cell
        remove_frommesh!(mesh, edge)
        remove_frommesh!(mesh, edge.twinEdge)
        if keepvert
            keepvert = false
            continue
        end
        remove_frommesh!(mesh, edge.originVertex)
    end
    remove_frommesh!(mesh, cell)
end


end # module
