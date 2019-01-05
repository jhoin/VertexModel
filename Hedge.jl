# This module holds the half edge information
module Hedge_class

export Hedge, newedgelen!, setedgeborder!, fillnexthedge

include("./Vertex.jl")
using .Vertex_class

# This object holds the half edge information
mutable struct Hedge <: AbstractDcel
    twinEdge::Hedge
    nextEdge::Hedge
    prevEdge::Hedge
    originVertex
    containCell
    edgeLen::Float64
    border::Bool

    Hedge() = new()
end

# Acessory function, it feills some of the fields in the next Hedge object
# Arguments: current edge, vertex object and cell object
# Return: DCEL object
function fillnexthedge(this, vert, cell)
    next = Hedge()
    next.containCell = cell
    next.prevEdge = this
    next.originVertex = vert
    return next
end # fillnexthedge

# Acessory function, it feills some of the fields in the first Hedge object
# Arguments: current edge, vertex object and cell object
# Return: DCEL object
function fillnexthedge(vert, cell)
    next = Hedge()
    next.containCell = cell
    next.originVertex = vert
    return next
end # fillnexthedge

# Find the edges located at the border in list of edges
# NOTE: Run only after the instantiation!!!!
# Arguments: list of edges
# Return: list of edge objects
function setedgeborder!(edge)
    edge.border = !isdefined(edge, :twinEdge)
end # assignborders

# Update edge length
# Arguments: edge object
# Return: edge object with length updated
function newedgelen!(edge)
    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex
    edge.edgeLen = distvertices(p1,p2)
end # newedgelen

end
