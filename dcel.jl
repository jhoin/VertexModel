# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future

module DCEL

abstract type AbstractDcel end

export importfromfile, newedgelen!, updatecell!, Dcel

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

# Create list of vertices
# Arguments: array of lines containing vertex coordinates
# Return: list of vertices objects
function getvertexlist(lines::Array{String,1},nPoints::Integer)
    pointCoords = Array{Float64}(undef,2)
    vertices = Array{Vertex}(undef,nPoints)
    for i in 1:nPoints
        pointCoords[:] = [parse(Float64,str) for str in split(popfirst!(lines))]
        vertex = Vertex()
        vertex.x = pointCoords[1]
        vertex.y = pointCoords[2]
        vertices[i] = vertex
    end # loop points

    return vertices
end # getvertexlist

# Find the twins of edges in list of edges
# Arguments: array of vertex index in each edge, number or edges
# Return: list of edge objects
function findtwinedges(provEdges::Vector{Hedge}, vertsInEdge::Array{Int64,2},nEdges::Integer)
    for i in 1:nEdges
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

    return provEdges
end # findtwinedges

# Acessory function, it feills some of the fields in the Hedge object
# Arguments:
# Return: DCEL object
function fillnexthedge(this::Hedge, vert::Vertex, cell::Cell)
    next = Hedge()
    next.containCell = cell
    next.prevEdge = this
    next.originVertex = vert
    return next
end # fillnexthedge

# Gets all verts in edge, the vert reference is the index on the vertex list
# Arguments:
# Return: Vector with the vertex in each edge
function getedgeverts(lines::Array{String,1}, info::Array{Int,1})

    vertsInCells = Vector{Int}(undef,info[2])
    vertsInEdge = Array{Int}(undef,info[1]*info[2],2)
    fill!(vertsInEdge,0)
    nEdges = 0
    lastVertexCell = 1
    for i in 1:info[1]
        vertsInCells[:] = [parse(Int,str) for str in split(lines[i])]

        # Get the edges starting with vertices in the middle
        for j in 1:info[2]-1
            if(vertsInCells[j+1] < 0 )
                lastVertexCell = j
                break
            end
            nEdges += 1
            vertsInEdge[nEdges,:] = [vertsInCells[j],vertsInCells[j+1]]
            println(vertsInEdge[nEdges,:])
        end # loop cell vertices

        # Get the edge starting at the last vertex of the cell
        nEdges += 1
        vertsInEdge[nEdges,:] = [vertsInCells[lastVertexCell],vertsInCells[1]]
        println(vertsInEdge[nEdges,:])
    end # loop cell list

    return vertsInEdge[1:nEdges,:]
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
    vertices = getvertexlist(allLines,nPoints) # allLines is modified inside function!!!

    # Put edges and cells inside arrays
    cellListInfo = [parse(Int,str) for str in split(popfirst!(allLines))]
    vertsInEdge = getedgeverts(allLines,cellListInfo)
    vertsInCells = Vector{Int}(undef,cellListInfo[2])

    edges = Array{Hedge}(undef,size(vertsInEdge,1))
    cells = Array{Cell}(undef,cellListInfo[1])
    nEdges = 1
    lastVertexCell = 1
    for i in 1:cellListInfo[1]
        vertsInCells[:] = [parse(Int,str) for str in split(popfirst!(allLines))]

        # Get the edge starting at the first vertex of the cell
        thisCell = Cell()
        firstEdge = Hedge()
        firstEdge.originVertex = vertices[vertsInCells[1]]
        vertices[vertsInCells[1]].leavingEdge = firstEdge
        firstEdge.containCell = thisCell
        thisCell.incEdge = firstEdge
        thisEdge = firstEdge
        cells[i] = thisCell

        # Get the edges starting with vertices in the middle
        for j in 2:cellListInfo[2]-1
            if(vertsInCells[j+1] < 0 )
                lastVertexCell = j
                break
            end
            nextEdge = fillnexthedge(thisEdge, vertices[vertsInCells[j]], thisCell)
            thisEdge.nextEdge = nextEdge

            vertices[vertsInCells[j]].leavingEdge = nextEdge
            nEdges += 1
            edges[nEdges] = nextEdge
            thisEdge = nextEdge
        end # loop cell vertices,

        # Get the edge starting at the last vertex of the cell
        lastEdge = fillnexthedge(thisEdge, vertices[vertsInCells[lastVertexCell]], thisCell)
        lastEdge.nextEdge = firstEdge
        firstEdge.prevEdge = lastEdge
        thisEdge.nextEdge = lastEdge

        vertices[vertsInCells[lastVertexCell]].leavingEdge = lastEdge

        nEdges += 1
        edges[nEdges] = lastEdge
        edges[1] = firstEdge
    end # loop cell list

    # Set the twin edges
    findtwinedges(edges,vertsInEdge,nEdges)

    return Dcel(vertices,edges,cells)
end # importfromfile

# Update edge length
# Arguments: edge object
# Return: edge object with length updated
# TODO: move to an edge module
function newedgelen!(edge::Hedge)
    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex
    edge.edgeLen = distvertices(p1,p2)
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

# Update cell measures
# Arguments: a cell object
# Return: float number storing the distance
# TODO: move to a cell module
function updatecell!(cell::Cell)
    edge = cell.incEdge
    p1 = edge.originVertex
    sumArea = 0.0
    sumPerim = 0.0
    first = p1
    p2 = edge.nextEdge.originVertex

    # loop over vertex and get the measurements
    while true
        sumArea += abs(p1.x*p2.y - p2.x*p1.y)
        sumPerim += distvertices(p1,p2)
        p1 = p2
        edge = edge.nextEdge
        if(p2 == first) break end
        p2 = edge.nextEdge.originVertex
    end

    # Update the values
    cell.areaCell = sumArea
    cell.perimCell = sumPerim

    return cell
end # updatecell

end # module

using .DCEL

system = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex2.points")

# Update the measurements
for i in 1:length(system.listEdge)
    newedgelen!(system.listEdge[i])
end

for i in 1:length(system.listCell)
    updatecell!(system.listCell[i])
end

#print(system.listEdge[1].twinEdge.edgeLen)

# TODO:
# Plotting operations for mesh
# FIX MISSING INFO BUG!!!!
# Export mesh to file
