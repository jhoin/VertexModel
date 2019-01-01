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

export importfromfile, Dcel

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

    #Hedge{T}() where {T} = new{T}()


    #Hedge(originVertex) = new{T}(originVertex,undef,undef,undef,undef,0.0)
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
        firstEdge.edgeLen = 1.0
        thisCell = Cell()
        thisCell.incEdge = firstEdge
        firstEdge.containCell = thisCell
        thisEdge = firstEdge


        cells[i] = thisCell
        vertsInEdge[1,:] = vertsInCells[1:2]
        #vertsInEdge[nEdges,2] = vertsInCells[2]

        # Get the edges starting with vertices in the middle
        for j in 2:cellListInfo[2]-1
            if(vertsInCells[j+1] < 0)
                lastVertexCell = j
                break
            end

            nextEdge = Hedge()
            nextEdge.originVertex = vertices[vertsInCells[j]]
            nextEdge.containCell = thisCell
            nextEdge.prevEdge = thisEdge
            nextEdge.edgeLen = 2.0
            thisEdge.nextEdge = nextEdge
            thisEdge = nextEdge
            nEdges += 1
            provEdges[nEdges] = thisEdge
            vertsInEdge[nEdges,:] = [vertsInCells[j],vertsInCells[j+1]]
        end # loop cell vertices

        # Get the edge starting at the last vertex of the cell
        lastEdge = Hedge()
        lastEdge.originVertex = vertices[vertsInCells[lastVertexCell]]
        lastEdge.containCell = thisCell
        lastEdge.edgeLen = 3.0
        lastEdge.prevEdge = thisEdge
        lastEdge.nextEdge = firstEdge
        firstEdge.prevEdge = lastEdge

        nEdges += 1
        provEdges[nEdges] = lastEdge
        provEdges[1] = firstEdge
        vertsInEdge[nEdges,:] = [vertsInCells[lastVertexCell],vertsInCells[1]]
    end # loop cell list

    # Set the twin edges
    for i in 1:nEdges
        verts = reverse(vertsInEdge[i,:])
        for j in 1:nEdges
            all(verts == 0) && break
            i == j && continue
            found = all(vertsInEdge[j,:] == verts)
            #println(verts,"--->",vertsInEdge[j,:])
            if(found)
                provEdges[i].twinEdge = provEdges[j]
                break
            end
        end
    end # loop edges

    return Dcel(vertices,provEdges,cells)
end # loop function

end # module

using .DCEL

system = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex2.points")
for i in 1:length(system.listEdge)
    edge = system.listEdge[i]
    println(edge.edgeLen)
end

# TODO:
# Import mesh from file
# Update mesh measurements (area, perimeter, lenght...)
# Export mesh to file
# Plotting operations for mesh
