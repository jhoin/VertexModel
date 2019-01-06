# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future

module DCEL

include("Cell.jl")
include("Vertex.jl")
include("Hedge.jl")

using .Cell_class
using .Vertex_class
using .Hedge_class

export importfromfile, updatesystem!, t1transition!, Dcel

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
        vertex = createvertex(pointCoords)
        vertex.leavingEdges = Vector{Hedge}(undef,3)
        vertices[i] = vertex
    end # loop points

    return vertices
end # getvertexlist

# Find the twins of edges in list of edges
# Arguments: array of vertex index in each edge, number or edges
# Return: list of edge objects
function findtwinedges!(provEdges::Vector{Hedge}, vertsInEdge::Array{Int64,2},nEdges::Integer)
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
end # findtwinedges

# Gets all verts in edge, the vert reference is the index on the vertex list
# Arguments: connectivity matrix (from file) and array with # of cells and verts in cells
# Return: Vector with the vertex in each edge
function getedgeverts(lines::Array{String,1}, info::Array{Int,1})

    vertsInCells = Vector{Int}(undef,info[2])
    vertsInEdge = Array{Int}(undef,info[1]*info[2],2)
    fill!(vertsInEdge,0)
    nEdges = 1
    lastVertexCell = 1
    for i in 1:info[1]
        vertsInCells[:] = [parse(Int,str) for str in split(lines[i])]

        # Get the edges starting with vertices in the middle
        for j in 1:info[2]-1
            if(vertsInCells[j+1] < 0 )
                lastVertexCell = j
                break
            end
            vertsInEdge[nEdges,:] = [vertsInCells[j],vertsInCells[j+1]]
            #println(vertsInEdge[nEdges,:])
            nEdges += 1
        end # loop cell vertices

        # Get the edge starting at the last vertex of the cell
        vertsInEdge[nEdges,:] = [vertsInCells[lastVertexCell],vertsInCells[1]]
        #println(vertsInEdge[nEdges,:])
    end # loop cell list

    return vertsInEdge[1:nEdges,:]
end # getedgeverts

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
    nEdges,lastVertexCell = 1,1
    for i in 1:cellListInfo[1]
        vertsInCells[:] = [parse(Int,str) for str in split(popfirst!(allLines))]

        # Get the edge starting at the first vertex of the cell
        thisCell = Cell()
        firstEdge = fillnexthedge(vertices[vertsInCells[1]], thisCell)
        setleavingedge!(vertices[vertsInCells[1]],firstEdge)
        firstEdgeIndex = nEdges

        thisCell.incEdge = firstEdge
        thisEdge = firstEdge
        cells[i] = thisCell

        # Get the edges starting with vertices in the middle
        cellListInfo,vertsInCells,edges,nEdges
        for j in 2:cellListInfo[2]
            if(j == cellListInfo[2] || vertsInCells[j+1] < 0)
                lastVertexCell = j
                break
            end
            nextEdge = fillnexthedge(thisEdge, vertices[vertsInCells[j]], thisCell)
            thisEdge.nextEdge = nextEdge
            setleavingedge!(vertices[vertsInCells[j]],nextEdge)
            nEdges += 1
            edges[nEdges] = nextEdge
            thisEdge = nextEdge
        end # loop cell vertices,

        # Get the edge starting at the last vertex of the cell
        lastEdge = fillnexthedge(thisEdge, vertices[vertsInCells[lastVertexCell]], thisCell)
        lastEdge.nextEdge = firstEdge
        firstEdge.prevEdge = lastEdge
        thisEdge.nextEdge = lastEdge
        setleavingedge!(vertices[vertsInCells[lastVertexCell]],lastEdge)
        nEdges += 1
        edges[nEdges] = lastEdge
        edges[firstEdgeIndex] = firstEdge
    end # loop cell list

    # Set the twin edges
    findtwinedges!(edges,vertsInEdge,nEdges)

    return Dcel(vertices,edges,cells)
end # importfromfile

# Update the measures of the mesh, namely cell areas and perimeters and edge lens
# Arguments: DCEL object
# Return: DCEL object
function updatesystem!(system)

    # Update cells
    for i in 1:length(system.listCell)
        updatecell!(system.listCell[i])
        if(system.listCell[i].areaCell < 0.0)
            invertcell!(system.listCell[i])
            updatecell!(system.listCell[i])
        end
    end

    # Update edges
    for i in 1:length(system.listEdge)
        newedgelen!(system.listEdge[i])
        #println(system.listEdge[i].edgeLen)
        setedgeborder!(system.listEdge[i])
        println(system.listEdge[i].border)
    end

    # Look for topological changes
    #t1transition!(system.listEdge[149])
    # for i in 1:length(system.listEdge)
    #     system.listEdge[i].border && continue
    #     if system.listEdge[i].edgeLen < 0.1
    #         t1transition!(system.listEdge[i])
    #     end
    # end
end # updatesystem

# Perform t1 transition
# Arguments:edge object
# Return: edge object updated
function t1transition!(focal_edge)

    # Perform the topological changes for the focal edge
    t1topology!(focal_edge)

    # Perform the topological changes for the focal twin
    t1topology!(focal_edge.twinEdge)

    # Rotate points
    p1 = focal_edge.originVertex
    p2 = focal_edge.nextEdge.originVertex

end # t1transition

# Make the topological changes for the t1 transition
# Arguments:edge object
# Return: edge object updated
function t1topology!(focal_edge)

    # Check if edge is not incident edge of cell
    focal_cell = focal_edge.containCell
    if focal_cell.incEdge == focal_edge
        focal_edge.incEdge = focal_edge.nextEdge
    end

    # Edge operations
    edge = focal_edge.prevEdge
    println(edge.border)
    edgeremove!(focal_edge)
    addedgeat!(focal_edge, edge)
end # t1topolgy

end # module

using .DCEL

system = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex2.points")
updatesystem!(system)

# TODO:
# Split the file into modules
# Plotting operations for mesh
# Export mesh to file
