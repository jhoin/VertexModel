# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future

# This object holds the cell information
mutable struct Cell{T}
    incEdge::T
    perimCell::Float64
    areaCell::Float64
    centroidCell::Vector{Float64}
end

# This object holds the half edge information
mutable struct Hedge{T}
    originVertex::T
    twinEdge::Hedge
    nextEdge::Hedge
    prevEdge::Hedge
    containCell::Cell
    edgeLen::Float64
end

# This object hold the vertex information
mutable struct Vertex{T}
    x::T
    y::T
    leavingEdge::Hedge
end

# contains the whole mesh
mutable struct Dcel
    listVert::Vector{Vertex}
    listEdge::Vector{Hedge}
    listCell::Vector{Cell}
end
