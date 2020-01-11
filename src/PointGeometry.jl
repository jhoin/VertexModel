module PointGeometry

using ..Meshes

# Calculate the distance between two vertices
# Arguments: two Vertex objects
# Return: float number storing the distance
# TODO: move to a vertex module
function pointdistance(p1::Vertex, p2::Vertex)
    x1 = p1.x
    y1 = p1.y
    x2 = p2.x
    y2 = p2.y

    return sqrt((x2-x1)^2 + (y2-y1)^2)
end # distvertices

# Rotate an edge around a middle point
# Arguments: two vertices
# Return: vertices with positions updated
function rotatepoints!(p1::Vertex, p2::Vertex)

    # midpoint
    cx = (p1.x + p2.x) / 2.0
    cy = (p1.y + p2.y) / 2.0

    # Change vertex coordinates
    theta = -90.0* π/360.
    x1 = (  (p1.x - cx) * cos(theta) + (p1.y - cy) * sin(theta) ) + cx
    y1 = ( -(p1.x - cx) * sin(theta) + (p1.y - cy) * cos(theta) ) + cy

    x2 = (  (p2.x - cx) * cos(theta) + (p2.y - cy) * sin(theta) ) + cx
    y2 = ( -(p2.x - cx) * sin(theta) + (p2.y - cy) * cos(theta) ) + cy

    # update the objects
    p1.x,p1.y = x1,y1
    p2.x,p2.y = x2,y2
end # rotatepoints
export rotatepoints!

# Rotate an edge around the first point
# Arguments: two vertices
# Return: vertices with positions updated
function rotatepoint(p1::Vertex, p2::Vertex)

    xy1 = [p1.x, p1.y]
    ptranslated = [p2.x-p1.x, p2.y-p1.y]
    angle_matrix = [cos(π) -sin(π); sin(π) cos(π)]

    # Change vertex coordinates
    xy2 = angle_matrix*ptranslated + xy1

    cx = (p1.x + p2.x) / 2.0
    cy = (p1.y + p2.y) / 2.0
    theta = 1.5
    x2 = (  (p2.x - p1.x) * cos(theta) + (p2.y - p1.x) * sin(theta) ) + p1.x
    y2 = ( -(p2.x - p1.y) * sin(theta) + (p2.y - p1.y) * cos(theta) ) + p1.y
    moveAngle = atan((p1.y - p2.y) / (p1.x - p2.x))

    # update the object
    return xy2
end # rotatepoints
export rotatepoint

function intersection(p1::Vertex, p2::Vertex, p3::Vertex, p4::Vertex)
    a1 = p2.y - p1.y
    b1 = p1.x - p2.x
    c1 = a1 * p1.x + b1 * p1.y

    a2 = p4.y - p3.y
    b2 = p3.x - p4.x
    c2 = a2 * p3.x + b2 * p3.y

    delta = a1 * b2 - a2 * b1

    # If lines are parallel, intersection point will contain infinite values
    return Vertex([ (b2 * c1 - b1 * c2) / delta, (a1 * c2 - a2 * c1) / delta])
end # intersection
export intersection

function intersection(line1::Array{Float64,2}, line2::Array{Float64,2})
    a1 = line1[2,2] - line1[1,2]
    b1 = line1[1,1] - line1[2,1]
    c1 = a1 * line1[1,1] + b1 * line1[1,2]

    a2 = line2[2,2] - line2[1,2]
    b2 = line2[1,1] - line2[2,1]
    c2 = a2 * line2[1,1] + b2 * line2[1,2]

    delta = a1 * b2 - a2 * b1
    return [(b2 * c1 - b1 * c2) / delta, (a1 * c2 - a2 * c1) / delta]
end
export intersection

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
export is_between

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
    return Vertex([x, y])
end # lineintersect

function lineintersect(p1::Vertex, p2::Vertex, p3::Vertex, p4::Vertex)
    x_num = (p1.x*p2.y - p1.y*p2.x)*(p3.x - p4.x) - (p1.x - p2.x)*(p3.x*p4.y - p3.y*p4.x)
    y_num = (p1.x*p2.y - p1.y*p2.x)*(p3.y - p4.y) - (p1.y - p2.y)*(p3.x*p4.y - p3.y*p4.x)
    den = (p1.x - p2.x)*(p3.y - p4.y) - (p1.y - p2.y)*(p3.x - p4.x)
    return Vertex([x_num/den, y_num/den])
end

# Check if point lies in an edge
# NOTE: it is assumed the point lies in the line equation of the edge
# Arguments: edge object and point
# Return: true or false
function isvert_inedge(edge::Hedge, vert::Vertex, min::Float64)
    distAC = pointdistance(edge.originVertex, vert)
    distBC = pointdistance(edge.nextEdge.originVertex, vert)
    return (edge.edgeLen - distAC + distBC ) > min
end # isvert_inedge

# Extend an edge and its twin in 50% the minimum edge lenght
# Arguments: edge
# Return: edge and twin updated
function extend_edge(edge::Hedge, minlen)

    # Get points
    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex

    # Move the points in the direction given by the above angle
    moveAngle = atan((p1.y - p2.y) / (p1.x - p2.x))
    p2.x = p2.x + 0.3*minlen*cos(moveAngle)
    p2.y = p2.y + 0.3*minlen*sin(moveAngle)
    p1.x = p1.x + 0.3*minlen*cos(moveAngle)
    p1.y = p1.y + 0.3*minlen*sin(moveAngle)

    # Update edge lenght
    newedgelen!(edge)
    newedgelen!(edge.twinEdge)
end # extend_edge
export extend_edge

# Shorten an edge by 30%
# Arguments: eedge
# Return: edge and twin updated
function shorten_edge!(edge::Hedge)
    # Get points
    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex

    # Move the points in the direction given by the above angle
    moveAngle = atan((p1.y - p2.y) / (p1.x - p2.x))
    p2.x = p2.x - 0.05*cos(moveAngle)
    p2.y = p2.y - 0.05*sin(moveAngle)

    # Update edge lenght
    newedgelen!(edge)
    newedgelen!(edge.twinEdge)
end # shorten_edge!

# Calculate the shortest axis and return a point on that axis
# Arguments: cell
# Return: vertex on the shortest axis
function get_shortaxis(cell::Cell)

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
    axis_angle = inertia_eigen.vectors[idx,:]

    #centroid_displace = [axis_angle[1]+cell.centroid.x, axis_angle[2]+cell.centroid.y]
    centroid_displace = [axis_angle[1]+1.0, axis_angle[2]+1.0]
    return Vertex(centroid_displace)
end # get_shortaxis

# Get short axis as the perpendicar line to the farthest point
# Arguments: cell object
# Return: Vertex in the direction of the shortest axis
function shortest_axis(cell::Cell)
    farthest = Vertex([0.0, 0.0])
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

    return Vertex([x,y])
end # shortest_axis
export shortest_axis

"""
    attempt_move!(vert::Vertex, radius::Float64)
Move a vertex in random angle at a radius distance.
"""
function attempt_move!(vert::Vertex, radius::Float64)
    rnd1 = rand()
    rnd2 = rand()

    # If one is bigger than the other, swap them
    if rnd2 > rnd1 rnd2, rnd1 = rnd1, rnd2 end

    # Moooooove the points
    angle = 360*rand()
    vert.x = vert.x + rnd2*radius*cos(angle*π*180.0)
    vert.y = vert.y + rnd2*radius*sin(angle*π*180.0)
    #vert.x = vert.x + rnd2*radius*cos(2.0*π*rnd1/rnd2)
    #vert.y = vert.y + rnd2*radius*sin(2.0*π*rnd1/rnd2)
end
export attempt_move!

end # end of module
