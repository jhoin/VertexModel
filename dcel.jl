# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future

module DCEL

include("Cell.jl")
include("Vertex.jl")
include("Hedge.jl")

using .Cell_class
using .Vertex_class
using .Hedge_class
using DelimitedFiles

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

# Gets all verts in edge, the vert reference is the index on the vertex list
# Arguments: connectivity matrix (from file) and array with # of cells and verts in cells
# Return: Vector with the vertex in each edge
function getedgeverts(lines::Array{String,1}, info::Array{Int,1})

    vertsInCells = Vector{Int}(undef,info[2])
    vertsInEdge = Array{Int}(undef,info[1]*info[2],2)
    fill!(vertsInEdge,0)
    nEdges,lastVertexCell = 1,1
    for i in 1:info[1]
        vertsInCells[:] = [parse(Int,str) for str in split(lines[i])]

        # Get the edges starting with vertices in the middle
        for j in 1:info[2]-1
            if(vertsInCells[j+1] < 0 )
                lastVertexCell = j
                break
            end
            vertsInEdge[nEdges,:] = [vertsInCells[j],vertsInCells[j+1]]
            nEdges += 1
        end # loop cell vertices

        # Get the edge starting at the last vertex of the cell
        vertsInEdge[nEdges,:] = [vertsInCells[lastVertexCell],vertsInCells[1]]
    end # loop cell list

    return vertsInEdge[1:nEdges,:]
end # getedgeverts

# Gets all verts in a cell, the vert reference is the index on the vertex list
# Arguments: connectivity matrix (from file)
# Return: Vector of vectors with the vertex in each cell for all cells
function getcellverts(lines, info)
    cellverts = []
    vertsInCells = Vector{Int}(undef,info[2])
    for i in 1:info[1]
        vertsInCells[:] = [parse(Int,str) for str in split(lines[i])]
        push!(cellverts,vertsInCells[vertsInCells .> 0])
    end
    return cellverts
end # getcellverts

# Gets the number of verts in each cell
# Arguments: connectivity matrix
# Return: vector with cell sizes
function getcellsizes(cellverts)
    cellsizes = Vector{Int64}(undef,size(cellverts))
    for i in 1:size(cellverts)[1]
        cellsizes[i] = size(cellverts[i])[1]
    end
    return cellsizes
end # getcellsizes

# Build the verts of edge table, and then
# Arguments: vertices in each cell
# Return: verts in each edge foe all edges
function getedgeverts(cellverts)
    edgeverts = []
    for i in 1:size(cellverts,1)
        for j in 1:size(cellverts[i],1)-1
            push!(edgeverts,[cellverts[i][j], cellverts[i][j+1]])
        end
        push!(edgeverts,[cellverts[i][size(cellverts[i],1)], cellverts[i][1]])
    end
    return edgeverts
end # getedgeverts

# Create list of edges and assign the origin vertex to all of them
# Arguments: edge list, vertex list and vertsInEdge table
# Return: edge list updated
function findorigin!(edges, vertices, vertsInEdge)
    for i in 1:size(vertsInEdge,1)
        edges[i] = Hedge()
        vert = vertices[vertsInEdge[i,1]]
        edges[i].originVertex = vertices[vertsInEdge[i,1]][1]
    end
end # findorigin

# Assign the next edge field in list of edges
# NOTE: This function requires a list of the number of vertices per cell
# Arguments: edge list, list of cell sizes and vertsInEdge table
# Return: edge list updated
function findnext!(edges, vertsInEdge, cellsizes)
    edgecell = 1 # count the number of edges that I intereact in a cell
    icell = 1 # count the cells
    for i in 1:size(vertsInEdge,1)
        if edgecell == cellsizes[icell]
            edges[i].nextEdge = edges[i-cellsizes[icell]+1]
            icell += 1
            edgecell = 1
            continue
        end
        edges[i].nextEdge = edges[i+1]
        edgecell += 1
    end
end # findnext

# Assign the previous edge field in list of edges
# NOTE: This function requires a list of the number of vertices per cell
# Arguments: edge list, list of cell sizes and vertsInEdge table
# Return: edge list updated
function findprevious!(edges, vertsInEdge, cellsizes)
    edgecell = 1 # count the number of edges that I intereact in a given cell
    icell = 1  # count the cells
    for i in reverse(1:size(vertsInEdge,1))
        if edgecell == cellsizes[icell]
            edges[i].prevEdge = edges[i+cellsizes[icell]-1]
            icell += 1
            edgecell = 1
            continue
        end
        edges[i].prevEdge = edges[i-1]
        edgecell += 1
    end
end # findnext

# Find the twins of edges in list of edges
# Arguments: array of vertex index in each edge, number or edges
# Return: list of edge objects
function findtwinedges!(provEdges, vertsInEdge,nEdges)
    ntwins = 0
    for i in 1:nEdges
        verts = reverse(vertsInEdge[i,:][1])
        for j in 1:nEdges
            all(verts == 0) && break
            i == j && continue
            found = all(vertsInEdge[j,:][1] == verts)
            if(found)
                provEdges[i].twinEdge = provEdges[j]
                ntwins += 1
                break
            end
        end
    end
end # findtwinedges!

# Find the container cell of each edge
# Arguments: edges list, cell list and cell sizes
# Return: list of edge objects updated
function findcontainercell!(edges, cells, cellsizes)
    edgeid = 1
    cellid = 1
    for i in 1:size(edges,1)
        if edgeid == size(cellsizes[cellid],1)
            edgeid = 1
            cellid += 1
        end
        edges[edgeid].containCell = cells[cellid]
        edgeid += 1
    end
end

# Create the list of cells and assign the incident edge for each cell
# Arguments: cell list, edge list and list of cell sizes
# Return: cell list updated
function createlistcell(edges, cellsizes)
    cells =  Array{Cell}(undef,size(cellsizes,1))
    edgeid = 1
    for i in 1:size(cellsizes,1)
        cells[i] = Cell()
        cells[i].incEdge = edges[edgeid]
        edgeid += size(cellsizes[i],1)-1
    end
    return cells
end # createlistcell

# Import mesh from file, and represent it as a DCEL object
# Arguments: file path
# Return: DCEL object
function importfromfile(filePath::String)
    inFile = open(filePath)
    allLines = readlines(inFile)

    # Get the point coordinates
    nPoints = parse(Int,popfirst!(allLines))
    vertices = getvertexlist(allLines,nPoints) # allLines is modified inside function!!!

    # Put edges and cells inside arrays
    cellListInfo = [parse(Int,str) for str in split(popfirst!(allLines))]
    cellverts = getcellverts(allLines,cellListInfo)
    cellsizes = getcellsizes(cellverts)
    vertsInEdge = getedgeverts(cellverts)

    # Create list of edges
    edges = Array{Hedge}(undef,size(vertsInEdge,1))
    findorigin!(edges, vertices, vertsInEdge)
    findnext!(edges, vertsInEdge, cellsizes)
    findprevious!(edges, vertsInEdge, cellsizes)
    findtwinedges!(edges,vertsInEdge,size(vertsInEdge,1))

    # Create list of cells
    cells = createlistcell(edges,cellsizes)

    # connect cells and edges lists
    findcontainercell!(edges,cells,cellsizes)

    return Dcel(vertices,edges,cells)
end # importfromfile

# Update the measures of the mesh, namely cell areas and perimeters and edge lens
# Arguments: DCEL object
# Return: DCEL object
function updatesystem!(system::Dcel)

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
        setedgeborder!(system.listEdge[i])
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

# Export Dcel object to file
# Arguments: file path
# Return: DCEL object
function exporttofile(mesh::Dcel)

    # Create the coordinates list
    coords = Array{Float64}(undef, size(mesh.listVert,1), 2)
    for i in 1:size(mesh.listVert,1)
        coords[i,:] = [mesh.listVert[i].x, mesh.listVert[i].y]
    end

    # Create the connectivity list
    celltab = []
    for i in 1:size(mesh.listCell, 1)
        verts_cells = Cell_class.getcellverts(mesh.listCell[i])
        println(size(verts_cells,1))
        vertids = [findthisvert(verts_cells[vert], mesh.listVert) for vert in 1:size(verts_cells,1)]
        #println(size(verts_cells,1))
        append!(celltab, vertids)
    end

    #println(celltab)

    open("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex3_out.txt", "w") do out
        write(out, size(mesh.listVert,1))
        writedlm(out, coords)
        write(out, size(mesh.listCell,1))
        writedlm(out, celltab)
    end
end # exporttofile
export exporttofile

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
    edgeremove!(focal_edge)
    addedgeat!(focal_edge, edge)
end # t1topolgy

end # module

using .DCEL

system = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex3.points")
updatesystem!(system)
exporttofile(system)

# TODO:
# Plotting operations for mesh
# Export mesh to file
