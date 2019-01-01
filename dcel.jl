# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future

# include("Vertex.jl")
# include("Hedge.jl")
# include("Cell.jl")

#import .Vertex
#import .Hedge
#import .Cell

module DCEL

abstract type AbstractDcel end

export importfromfile, newedgelen, Dcel

# This object holds the cell information
mutable struct Cell <: AbstractDcel
    incEdge::AbstractDcel
    perimCell::Float64
    areaCell::Float64

    Cell() = new()
end


# This object holds the half edge information
mutable struct Hedge <: AbstractDcel
    originVertex::AbstractDcel
    twinEdge::Hedge
    nextEdge::Hedge
    prevEdge::Hedge
    containCell::Cell
    edgeLen::Float64

    Hedge() = new()
end

# This object hold the vertex information
mutable struct Vertex <: AbstractDcel
    x::Float64
    y::Float64
    leavingEdge::Hedge

    Vertex() = new()
end

# contains the whole mesh
mutable struct Dcel
    listVert::Vector{Vertex}
    listEdge::Vector{Hedge}
    listCell::Vector{Cell}
end

# Import mesh from file, and represent it as a DCEL object
# Arguments: file path
# Return: DCEL object
# TODO: break this function down... it's too long
function importfromfile(filePath::String)
    inFile = open(filePath)
    allLines = readlines(inFile)

    # Get the point coordinates
    nPoints = parse(Int,popfirst!(allLines))
    pointCoords = Array{Float64}(undef,2)
    vertices = Array{Vertex}(undef,nPoints)
    for i in 1:nPoints
        pointCoords[:] = [parse(Float64,str) for str in split(popfirst!(allLines))]
        vertex = Vertex()
        vertex.x = pointCoords[1]
        vertex.y = pointCoords[2]
        vertices[i] = vertex
    end # loop points

    # Put edges and cells inside arrays
    cellListInfo = [parse(Int,str) for str in split(popfirst!(allLines))]
    vertsInCells = Vector{Int}(undef,cellListInfo[2])
    vertsInEdge = Array{Int}(undef,cellListInfo[1]*cellListInfo[2],2)
    fill!(vertsInEdge,0)
    provEdges = Array{Hedge}(undef,cellListInfo[1]*cellListInfo[2])
    cells = Array{Cell}(undef,cellListInfo[1])
    nEdges = 1
    lastVertexCell = 1
    for i in 1:cellListInfo[1]
        vertsInCells[:] = [parse(Int,str) for str in split(popfirst!(allLines))]

        # Get the edge starting at the first vertex of the cell
        firstEdge = Hedge()
        firstEdge.originVertex = vertices[vertsInCells[1]]
        thisCell = Cell()
        thisCell.incEdge = firstEdge
        firstEdge.containCell = thisCell
        thisEdge = firstEdge

        cells[i] = thisCell
        vertsInEdge[1,:] = vertsInCells[1:2]

        # Get the edges starting with vertices in the middle
        for j in 2:cellListInfo[2]-1
            if(vertsInCells[j+1] < 0 )
                lastVertexCell = j
                break
            end

            nextEdge = Hedge()
            nextEdge.originVertex = vertices[vertsInCells[j]]
            nextEdge.containCell = thisCell
            nextEdge.prevEdge = thisEdge
            thisEdge.nextEdge = nextEdge
            nEdges += 1
            provEdges[nEdges] = nextEdge
            vertsInEdge[nEdges,:] = [vertsInCells[j],vertsInCells[j+1]]
            thisEdge = nextEdge
        end # loop cell vertices

        # Get the edge starting at the last vertex of the cell
        lastEdge = Hedge()
        lastEdge.originVertex = vertices[vertsInCells[lastVertexCell]]
        lastEdge.containCell = thisCell
        provEdges[nEdges].nextEdge = lastEdge
        lastEdge.prevEdge = provEdges[nEdges]
        lastEdge.nextEdge = firstEdge
        firstEdge.prevEdge = lastEdge

        nEdges += 1
        provEdges[nEdges] = lastEdge
        provEdges[1] = firstEdge
        vertsInEdge[nEdges,:] = [vertsInCells[lastVertexCell],vertsInCells[1]]
    end # loop cell list

    # Set the twin edges
    for i in 1:nEdges
        #println(vertsInEdge[i,:])
        verts = reverse(vertsInEdge[i,:])
        for j in 1:nEdges
            all(verts == 0) && break
            i == j && continue
            found = all(vertsInEdge[j,:] == verts)
            if(found)
                provEdges[i].twinEdge = provEdges[j]
                break
            end
        end
    end # loop edges

    # Clean edges
    edges = provEdges[1:nEdges]

    return Dcel(vertices,edges,cells)
end # importfromfile

# Update edge length
# Arguments: edge object
# Return: edge object with length updated
# TODO: move to an edge module
function newedgelen(edge::Hedge)

    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex
    edge.edgeLen = distvertices(p1,p2)

    #return edge
end # newedgelen

# Calculate the distance between two vertices
# Arguments: two Vertex objects
# Return: float number storing the distance
# TODO: move to a vertex module
function distvertices(p1::Vertex,p2::Vertex)
    x1 = p1.x
    y1 = p1.y

    x2 = p2.x
    y2 = p2.y

    return sqrt((x2-x1)^2 + (y2-y1)^2)
end # distvertices

end # module

using .DCEL

system = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex2.points")

# for i in 1:length(system.listVert)
#     vert = system.listVert[i]
#     println(vert.x," ",vert.y)
#
# end
# println("End of the vert coordinate")

#newedgelen(system.listEdge[135])
#print(size(system.listEdge))
for i in 1:length(system.listEdge)
    edge = system.listEdge[i]
    newedgelen(edge)
    println(system.listEdge[i].edgeLen)
end

# TODO:
# Import mesh from file
# Update mesh measurements (area, perimeter, lenght...)
# Export mesh to file
# Plotting operations for mesh
