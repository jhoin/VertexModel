# This module hold the vertex information
module Vertex_class
abstract type AbstractDcel end

export AbstractDcel, Vertex, createvertex,distvertices, setleavingedge!

mutable struct Vertex <: AbstractDcel
    x::Float64
    y::Float64
    leavingEdges

    Vertex() = new()
end

# Create a new Vertex object
# Arguments: coordinates
# Return: vertex object
function createvertex(pointCoords::Array{Float,1})
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
end
