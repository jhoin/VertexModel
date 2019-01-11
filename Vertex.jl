# This module hold the vertex information
module Vertex_class
abstract type AbstractDcel end

export AbstractDcel, Vertex, createvertex,distvertices, setleavingedge!, rotatepoints!

mutable struct Vertex <: AbstractDcel
    x::Float64
    y::Float64
    leavingEdges

    Vertex() = new()
end

x(vert::Vertex) = vert.x

# Create a new Vertex object
# Arguments: coordinates
# Return: vertex object
function createvertex(pointCoords)
    vert = Vertex()
    vert.x = pointCoords[1]
    vert.y = pointCoords[2]
    return vert
end

# Calculate the distance between two vertices
# Arguments: two Vertex objects
# Return: float number storing the distance
# TODO: move to a vertex module
function distvertices(p1,p2)
    x1 = p1.x
    y1 = p1.y
    x2 = p2.x
    y2 = p2.y

    return sqrt((x2-x1)^2 + (y2-y1)^2)
end # distvertices

# insert an edge to the list of edges leaving a vertex
# Arguments: vertex, hedge
# Return: updated vertex object
function setleavingedge!(vert,edge)
    for i in 1:3
        if !isassigned(vert.leavingEdges,i)
            vert.leavingEdges[i] = edge
        end
    end
end # setleavingedge!

# Find a vertex in a list of vertices and return its index
# Arguments: vertex, list of vertices
# Return: vert index or  if not found
function findthisvert(vert, vertices)
    found = 0 # index
    for i in 1:size(vertices,1)
        found = i
        vert == vertices[i] && break
        found = 0
    end
    return found
end # findthisvert
export findthisvert

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



end # module
