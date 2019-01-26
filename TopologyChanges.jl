module TopologyChanges

include("dcel.jl")
using .DCEL

# Check for topological changes in mesh
# Arguments: mesh
# Return: updated mesh
function update_topology!(mesh, minlen)

    # Check for t1 transitions
    for i in 1:size(mesh.listEdge,1)
        if mesh.listEdge[i].border
            continue
        end
        if mesh.listEdge[i].edgeLen < minlen
            t1transition!(mesh.listEdge[i], minlen)
            #break
        end
    end
end #update_topology
export update_topology!

# Perform t1 transition
# Arguments:edge object
# Return: edge object updated
function t1transition!(focal_edge,minlen)

    # Rotate points
    p1 = focal_edge.originVertex
    p2 = focal_edge.nextEdge.originVertex
    rotatepoints!(p1, p2)

    # Perform the topological changes for the focal edge
    t1topology!(focal_edge)

    # Perform the topological changes for the focal twin
    t1topology!(focal_edge.twinEdge)

    # Update edge length
    extend_edge(focal_edge, minlen)
end # t1transition

# Rotate an edge around a middle point
# Arguments: two vertices
# Return: vertices with positions updated
function rotatepoints!(p1,p2)

    # midpoint
    cx = (p1.x + p2.x) / 2.0
    cy = (p1.y + p2.y) / 2.0

    # Change vertex coordinates
    theta = -1.5
    x1 = (  (p1.x - cx) * cos(theta) + (p1.y - cy) * sin(theta) ) + cx
    y1 = ( -(p1.x - cx) * sin(theta) + (p1.y - cy) * cos(theta) ) + cy

    x2 = (  (p2.x - cx) * cos(theta) + (p2.y - cy) * sin(theta) ) + cx
    y2 = ( -(p2.x - cx) * sin(theta) + (p2.y - cy) * cos(theta) ) + cy

    # update the object
    p1.x,p1.y = x1,y1
    p2.x,p2.y = x2,y2
end # rotatepoints


# Extend an edge and its twin in 50% the minimum edge lenght
# Arguments: edge
# Return: edge and twin updated
function extend_edge(edge, minlen)

    # Get points
    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex

    # The projection of point(x1,x2) in a line paralle to x axis
    y_proj = p1.y
    x_proj = p2.x

    # Get the sides of a triangle formed by the
    moveAngle = atan((p1.y - p2.y) / (p1.x - p2.x))
    println("moveAngle: ", moveAngle)

    # Move the points in the direction given by the above angle
    p2.x = p2.x + 0.30*minlen*cos(moveAngle)
    p2.y = p2.y + 0.30*minlen*sin(moveAngle)

    # Update edge lenght
    newedgelen!(edge)
    println("Edge len: ",edge.edgeLen)
    newedgelen!(edge.twinEdge)
    println("Twin edge len: ",edge.twinEdge.edgeLen)
end # extend_edge

# Make the topological changes for the t1 transition
# Arguments:edge object
# Return: edge object updated
function t1topology!(focal_edge)
    println(focal_edge.edgeLen)

    # Check if edge is not incident edge of cell
    focal_cell = focal_edge.containCell
    if focal_cell.incEdge == focal_edge
        focal_cell.incEdge = focal_edge.nextEdge
    end

    # Edge operations
    edge = focal_edge.prevEdge.twinEdge
    edgeremove!(focal_edge, edge)
    addedgeat!(focal_edge, edge)
end # t1topolgy

# Remove edge from cell
# Arguments: edge object
# Return: edge object
function edgeremove!(edge, prevtwin)
    edge.prevEdge.nextEdge = edge.nextEdge
    edge.containCell = prevtwin.containCell
    #edge.nextEdge.originVertex = edge.twinEdge.originVertex
end # edgeremove

# Add edge into cell before indicated edge
# Arguments: edge to be added, and edge that will come after the given edge
# Return: edge object
function addedgeat!(focal_edge, next_edge)
    next_edge.originVertex = focal_edge.nextEdge.originVertex
    focal_edge.nextEdge = next_edge
    focal_edge.prevEdge = next_edge.prevEdge
    next_edge.prevEdge.nextEdge = focal_edge
    next_edge.prevEdge = focal_edge
end # addedgeat

# Perform a cell division
# Arguments: mesh object
# Return: mesh object updated
function cell_division!(mesh, cell)
    edge = cell.incEdge
    p1 = edge.originVertex
    first = p1
    p2 = edge.nextEdge.originVertex

    # Get the inertia tensor
    ixx,ixy, iyy = 0.0,0.0,0.0
    while true
        ixx += (p1.x*p2.y - p2.x*p1.y)*((p1.y)^2 + p1.y*p2.y + (p2.y)^2)
        iyy += (p1.x*p2.y - p2.x*p1.y)*((p1.x)^2 + p1.x*p2.x + (p2.x)^2)
        ixy += (p1.x*p2.y - p2.x*p1.y)*(p1.x*p2.y + 2.0*p1.x*p1.y + 2.0*p2.x*p2.y + p2.x*p1.y)
        p1 = p2
        edge = edge.nextEdge
        if(p2 == first) break end
        p2 = edge.nextEdge.originVertex
    end
    ixx = ixx/12.0
    iyy = iyy/12.0
    ixy = ixy/24.0
    inertia_matrix = [ixx -ixy; -ixy  iyy]
    display(inertia_matrix)
    #print("\n")
    axisdivision_angle = eigen(inertia_matrix)
    #println(axisdivision_angle.vectors[:,2])
    centroid_displace = [axisdivision_angle[1] + cell.cx, axisdivision_angle[1] + cell.cy]
    c2 = createvertex(centroid_displace)

end # cell_division!

# Write a vector to a file
# Arguments: vectorin the format: [xorig, yorig, ]
# Return: mesh object updated

end # module
