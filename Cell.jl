# This module holds the cell information

module Cell_class

include("./Vertex.jl")
using .Vertex_class
export Cell

# This object holds the cell information
mutable struct Cell <: AbstractDcel
    incEdge
    perimCell::Float64
    areaCell::Float64

    Cell() = new()
end

# Update cell measures
# Arguments: a cell object
# Return: float number storing the distance
function updatecell!(cell)
    edge = cell.incEdge
    p1 = edge.originVertex
    sumArea = 0.0
    sumPerim = 0.0
    first = p1
    p2 = edge.nextEdge.originVertex

    # loop over vertex and get the measurements
    while true
        sumArea += p1.x*p2.y - p2.x*p1.y
        sumPerim += distvertices(p1,p2)
        p1 = p2
        edge = edge.nextEdge
        if(p2 == first) break end
        p2 = edge.nextEdge.originVertex
    end

    # Update the values
    cell.areaCell = 0.5*sumArea
    cell.perimCell = sumPerim

    return cell
end # updatecell
export updatecell!

# Invert the order of the vertices in a cell
# Arguments: a cell object
# Return: DCEl object corrected
function invertcell!(cell)
    edge = cell.incEdge
    first = edge
    while true
        edge.nextEdge, edge.prevEdge = edge.prevEdge, edge.nextEdge
        edge = edge.prevEdge
        if(edge == first) break end
    end
end # invertcell
export invertcell!

# Transverse a cell and return the list of vertices
# Arguments: a cell object
# Return: list of vertices objects
function getcellverts(cell::Cell)
    verts = Vector{Vertex}()
    edge = cell.incEdge
    vert = edge.originVertex
    first = vert
    while true
        println(typeof(verts))
        push!(verts,vert)
        edge = edge.nextEdge
        vert = edge.nextEdge.originVertex
        if(vert == first) break end
    end
    return verts
end # getcellverts
export getcellverts


end
