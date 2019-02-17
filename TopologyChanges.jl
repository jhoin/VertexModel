module TopologyChanges

#include("dcel.jl")
using ..DCEL
using LinearAlgebra: eigen

# Check for topological changes in mesh
# Arguments: mesh
# Return: updated mesh
function update_topology!(mesh::Dcel, minlen::Float64)

    for i in 1:size(mesh.listCell,1)
       cell_division!(mesh, mesh.listCell[i])
    end

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
function rotatepoints!(p1::Vertex, p2::Vertex)

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
function extend_edge(edge::Hedge, minlen)

    # Get points
    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex

    # Move the points in the direction given by the above angle
    moveAngle = atan((p1.y - p2.y) / (p1.x - p2.x))
    p2.x = p2.x + 0.30*minlen*cos(moveAngle)
    p2.y = p2.y + 0.30*minlen*sin(moveAngle)

    # Update edge lenght
    newedgelen!(edge)
    newedgelen!(edge.twinEdge)
end # extend_edge

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
    edgeremove!(focal_edge, edge)
    addedgeat!(focal_edge, edge)
end # t1topolgy

# Remove edge from cell
# Arguments: edge object
# Return: edge object
function edgeremove!(edge::Hedge, prevtwin::Hedge)
    edge.prevEdge.nextEdge = edge.nextEdge
    edge.containCell = prevtwin.containCell
    #edge.nextEdge.originVertex = edge.twinEdge.originVertex
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

# Perform a cell division
# Arguments: mesh object
# Return: mesh object updated
function cell_division!(mesh::Dcel, cell::Cell)

    # Get the inertia tensor
    ixx,ixy, iyy = 0.0,0.0,0.0
    for edge in cell
        p1, p2 = edge.originVertex, edge.nextEdge.originVertex
        ixx += (p1.x*p2.y - p2.x*p1.y)*((p1.y)^2 + p1.y*p2.y + (p2.y)^2)
        iyy += (p1.x*p2.y - p2.x*p1.y)*((p1.x)^2 + p1.x*p2.x + (p2.x)^2)
        ixy += (p1.x*p2.y - p2.x*p1.y)*(p1.x*p2.y + 2.0*p1.x*p1.y + 2.0*p2.x*p2.y + p2.x*p1.y)
    end
    ixx = ixx/12.0
    iyy = iyy/12.0
    ixy = ixy/24.0
    inertia_matrix = [ixx -ixy; -ixy  iyy]
    inertia_eigen = eigen(inertia_matrix)
    idx = findmax(inertia_eigen.values)[2]
    println(inertia_eigen.vectors)
    axis_angle = inertia_eigen.vectors[idx,:]

    #centroid_displace = [axis_angle[1]+cell.centroid.x, axis_angle[2]+cell.centroid.y]
    centroid_displace = [axis_angle[1]+1.0, axis_angle[2]+1.0]
    println(centroid_displace)
    #
    return createvertex(centroid_displace)
end # get_shortaxis

# Get short axis as the perpendicar line to the farthest point
# Arguments: cell object
# Return: Vertex in the direction of the shortest axis
function shortest_axis(cell::Cell)
    farthest = Vertex()
    dist = 0.0
    for edge in cell
        p1 = edge.originVertex
        new_dist = distvertices(p1, cell.centroid)
        if dist < new_dist
            dist = new_dist
            farthest = p1
        end
    end
    angle_change = atan((farthest.x - cell.centroid.x)/(farthest.y - cell.centroid.y))
    x = cell.centroid.x + 0.1*cos(angle_change)
    y = cell.centroid.y + 0.1*sin(angle_change)

    return createvertex([x,y])
end # shortest_axis

#
function intersection(p1::Vertex, p2::Vertex, p3::Vertex, p4::Vertex)
    a1 = p2.y - p1.y
    b1 = p1.x - p2.x
    c1 = a1 * p1.x + b1 * p1.y

    a2 = p4.y - p3.y
    b2 = p3.x - p4.x
    c2 = a2 * p3.x + b2 * p3.y

    delta = a1 * b2 - a2 * b1

    # If lines are parallel, intersection point will contain infinite values
    return createvertex([ (b2 * c1 - b1 * c2) / delta, (a1 * c2 - a2 * c1) / delta])
end # intersection

# Check if point c lies on the line segment formed by a and b
# Arguments: vertex objects of each point
# Return: true (point lies in line) or false (does not lie)
function is_between(a::Vertex, b::Vertex, c::Vertex)
    crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y)

    # compare versus epsilon for floating point values, or != 0 if using integers
    if abs(crossproduct) > 0.005 return false end

    dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y)
    if dotproduct < 0 return false end

    squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y)
    if dotproduct > squaredlengthba return false end

    return true
end # is_between

# Compute slope of a line equation formed by two points
# Arguments: two vertex objects
# Return: a float containing the slope
function lineslope(p1::Vertex, p2::Vertex)
    return (p2.y - p1.y)/(p2.x - p1.x)
end # lineslope

# Compute intersect of the line equation in the y axis given a point and slope
# Arguments: vertex object and line slope
# Return: a float containing the intersect
function yintersect(p::Vertex, slope)
    return p.y - slope*p.x
end # yintersect

function yintersect(p1::Vertex, p2::Vertex)
    return (p1.x*p2.y - p1.y*p2.x) / (p1.x - p2.x)
end

# Compute the intersection point between two lines
# Arguments: slopes and intercepts of the lines
# Return: a vertex object of the point
function lineintersect(slope1, slope2, sect1, sect2)
    x = (sect2 - sect1) / (slope1 - slope2)
    y = (slope1*sect2 - slope2*sect1) / (slope1 - slope2)
    return createvertex([x, y])
end # lineintersect

function lineintersect(p1::Vertex, p2::Vertex, p3::Vertex, p4::Vertex)
    x_num = (p1.x*p2.y - p1.y*p2.x)*(p3.x - p4.x) - (p1.x - p2.x)*(p3.x*p4.y - p3.y*p4.x)
    y_num = (p1.x*p2.y - p1.y*p2.x)*(p3.y - p4.y) - (p1.y - p2.y)*(p3.x*p4.y - p3.y*p4.x)
    den = (p1.x - p2.x)*(p3.y - p4.y) - (p1.y - p2.y)*(p3.x - p4.x)
    return createvertex([x_num/den, y_num/den])
end

# Check if point lies in an edge
# NOTE: it is assumed the point lies in the line equation of the edge
# Arguments: edge object and point
# Return: true or false
function isvert_inedge(edge::Hedge, vert::Vertex, min::Float64)
    distAC = distvertices(edge.originVertex, vert)
    distBC = distvertices(edge.nextEdge.originVertex, vert)
    return (edge.edgeLen - distAC + distBC ) > min
end # isvert_inedge

end # module
