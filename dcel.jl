# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future

module DCEL

using DelimitedFiles

export importfromfile, updatesystem!, t1transition!, Dcel
abstract type AbstractDcel end


mutable struct Vertex <: AbstractDcel
    x::Float64
    y::Float64
    leavingEdges

    Vertex() = new()
end

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

# This object holds the cell information
mutable struct Cell <: AbstractDcel
    incEdge
    perimCell::Float64
    areaCell::Float64

    Cell() = new()
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
        vertex = createvertex(pointCoords)
        vertex.leavingEdges = Vector{Hedge}(undef,3)
        vertices[i] = vertex
    end # loop points

    return vertices
end # getvertexlist

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
export setleavingedge!

# Find a vertex in a list of vertices and return its index
# Arguments: vertex, list of vertices
# Return: vert index or  if not found
function findthisvert(vert, vertices)
    found = 0 # inde
    for i in 1:size(vertices,1)
        #found = i
        #println("Found: ",found)
        if vert == vertices[i]
            found = i
            #println("Found: ",found)
            break
        end
        #found = 0
    end
    return found
end # findthisvert

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

# Transverse a cell and return the list of vertices
# Arguments: a cell object
# Return: list of vertices objects
function getcellverts(cell::Cell)
    verts = Vector{Vertex}()
    edge = cell.incEdge
    #println("esge blongs to this cell: ", edge.containCell == cell)
    vert = edge.originVertex
    first = edge
    while true
        push!(verts,vert)
        edge = edge.nextEdge
        #println("edge len: ", edge.edgeLen)
        if(edge == first)
            #println("first edge len conditional: ", edge.edgeLen)
            break
        end
        vert = edge.originVertex
    end
    #println("n of verts cell: ",size(verts,1))
    return verts
end # getcellverts

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
        edges[i].containCell = cells[cellid]
        edgeid += 1
    end
end

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

# Remove edge from cell
# Arguments: edge object
# Return: edge object
function edgeremove!(edge)
    edge.prevEdge.nextEdge = edge.nextEdge
end # edgeremove

# Add edge into cell before indicated edge
# Arguments: edge to be added, and edge that will come after the given edge
# Return: edge object
function addedgeat!(focal_edge,next_edge)
    next_edge.originVertex = focal_edge.nextEdge.originVertex
    focal_edge.nextEdge = next_edge
    focal_edge.prevEdge = next_edge.prevEdge
    next_edge.prevEdge = focal_edge
end # addedgeat

# Create the list of cells and assign the incident edge for each cell
# Arguments: cell list, edge list and list of cell sizes
# Return: cell list updated
function createlistcell(edges, cellsizes)
    cells =  Array{Cell}(undef,size(cellsizes,1))
    edgeid = 1
    for i in 1:size(cellsizes,1)
        #println("index: ",i)
        cells[i] = Cell()
        cells[i].incEdge = edges[edgeid]
        edgeid += cellsizes[i]
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
        # if(system.listCell[i].areaCell < 0.0)
        #     invertcell!(system.listCell[i])
        #     updatecell!(system.listCell[i])
        # end
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
    for cell in 1:size(mesh.listCell,1)
        #println(mesh.listCell[cell] == mesh.listCell[cell+1])
        #println("Area: ",mesh.listCell[cell].areaCell)
        verts_cells = getcellverts(mesh.listCell[cell])
        vertids = [findthisvert(verts_cells[vert], mesh.listVert) for vert in 1:size(verts_cells,1)]
        #println("index: ",cell)
        push!(celltab, vertids)
    end

    #println(celltab)

    open("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex3_out.txt", "w") do out
        #println(size(mesh.listVert,1))
        write(out, size(mesh.listVert,1))
    open("/home/jhon/Documents/Projects/vertexModelJulia/results/ex3_out.txt", "w") do out
        nvert = size(mesh.listVert,1)
        write(out, "$nvert\n")
        writedlm(out, coords)
        ncell = size(mesh.listCell,1)
        write(out, "$ncell\n")
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

system = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/results/ex3.points")
updatesystem!(system)
exporttofile(system)

# TODO:
# Plotting operations for mesh
# Export mesh to file
