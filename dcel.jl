# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future
module Meshes

using DelimitedFiles
import Base: iterate, length, show
import Base:(==)

export importfromfile, updatesystem!, t1transition!
abstract type AbstractBoundary end
abstract type AbstractLine end
abstract type AbstractFace end

mutable struct Vertex
    x::Float64
    y::Float64
    id::Int64
    leavingEdges
    treatBoundary::AbstractBoundary

    Vertex(coords::Vector{Float64}) = createvertex(new(), coords)
end
export Vertex

# Boundary condition for vertices tied to springs
mutable struct Spring <: AbstractBoundary
    x::Float64
    y::Float64
    spring_const::Float64
    len::Float64
    vert::Vertex
    Spring(x, y, spring_const, len, vert) = new(x, y, spring_const, len, vert)
end
export Spring

# Boundary condition for immobile vertices
struct Immobile <: AbstractBoundary
    vert::Vertex
    Immobile(vert) = new(vert)
end
export Immobile

# Boundary condition for vertices without treatment, they can move
struct Mobile <: AbstractBoundary
    vert::Vertex
    Mobile(vert) = new(vert)
end
export Mobile

# This object holds the half edge information
mutable struct Hedge <: AbstractLine
    twinEdge::Hedge
    nextEdge::Hedge
    prevEdge::Hedge
    originVertex::Vertex
    containCell
    edgeLen::Float64
    border::Bool
    eqbond::Float64
    id::Int64
    Hedge() = create_hedge(new())
end
export Hedge

# This object holds the cell information
mutable struct Cell <: AbstractFace
    incEdge::Hedge
    perimCell::Float64
    areaCell::Float64
    celltype::Int64
    area_elast::Float64
    perim_elast::Float64
    centroid::Vertex
    id::Int64

    Cell() = new()
end
export Cell

# contains the whole mesh
struct Mesh
    vertices::Vector{Vertex}
    edges::Vector{Hedge}
    cells::Vector{Cell}
    lastid::Vector{Int64}
end
export Mesh

struct SubMesh
    vert::Vertex
    edges::Vector{Hedge}
    cells::Vector{Cell}
    unitlen::Float64
    #SubMesh(vert::Vertex, unitlen::Float64) = get_submesh(new(), vert, unitlen)
end
export SubMesh

function get_submesh(vert::Vertex, unitlen::Float64)
    edges = collect(skipmissing(vert.leavingEdges))
    cells = Vector{Cell}(undef, length(edges))
    for i in 1:length(edges)
        cells[i] = edges[i].containCell
    end
    submesh = SubMesh(vert, edges, cells, unitlen)
end
export get_submesh

(==)(x::Hedge, y::Hedge) = x.id == y.id
(==)(x::Cell, y::Cell) = x.id == y.id
(==)(x::Vertex, y::Vertex) = x.id == y.id

# Iterate over the vertices of a cell
# Arguments: array of lines containing vertex coordinates
# Return: list of vertices objects
function Base.iterate(iter::Cell)
    state = (iter.incEdge, false)
    return iterate(iter, state)
end

# Iterate over the vertices of a cell
# Arguments: array of lines containing vertex coordinates
# Return: list of vertices objects
function Base.iterate(iter::Cell, state)
    edge, pass_first = state

    if edge == iter.incEdge && pass_first
        return nothing
    end
    return edge, (edge.nextEdge, true)
end
Base.eltype(iter::Cell) = Hedge

function Base.length(cell::Cell)
    len = 0
    for edge in cell
        len += 1
    end
    return len
end

function Base.length(mesh::SubMesh)
    return length(mesh.edges)
end

function Base.show(mesh::Mesh)
    println("Number of vertices: ",length(mesh.vertices))
    println("Number of edges: ",length(mesh.edges))
    println("Number of cells: ",length(mesh.cells))
end

"""
    getvertexlist(lines::Array{String,1},nPoints::Integer)
Create list of vertices from a string taken from file
"""
function getvertexlist(lines::Array{String,1},nPoints::Integer)
    pointCoords = Array{Float64}(undef,2)
    vertices = Array{Vertex}(undef,nPoints)
    for i in 1:nPoints
        pointCoords[:] = [parse(Float64,str) for str in split(popfirst!(lines))]
        vertex = Vertex(pointCoords)
        vertex.id = i
        vertices[i] = vertex
    end # loop points
    return vertices
end # getvertexlist

"""
    getvertexlist(lines::Array{String,1}, nPoints::Integer)
Create list of vertices from coordinate matrix
"""
function getvertexlist(coords::Array{Float64,2})
    vertices = Array{Vertex}(undef,size(coords, 1))
    for i in 1:size(coords, 1)
        vertex = Vertex(coords[i,:])
        vertex.id = i
        vertices[i] = vertex
    end # loop points
    return vertices
end # getvertexlist

"""
    createvertex(pointCoords)
Create a new Vertex object from two coordinates and return it.
"""
function createvertex(vert::Vertex, pointCoords::Vector{Float64})
    vert.x = pointCoords[1]
    vert.y = pointCoords[2]
    vert.leavingEdges = Vector{Union{Missing,Hedge}}(missing,3)
    return vert
end
export createvertex

function create_hedge(edge::Hedge)
    #edge.originVertex = vert
    edge.containCell = Cell()
    return edge
end

# Calculate the distance between two vertices
# Arguments: two Vertex objects
# Return: float number storing the distance
# TODO: move to a vertex module
function distvertices(p1::Vertex, p2::Vertex)
    x1 = p1.x
    y1 = p1.y
    x2 = p2.x
    y2 = p2.y

    return sqrt((x2-x1)^2 + (y2-y1)^2)
end # distvertices
export distvertices

function update_spring!(spring::Spring)
    len = sqrt((spring.x - spring.x)^2 + (spring.y - spring.y)^2)
    spring.len = len
end
export update_spring!

"""
    setleavingedge!(vert::Vertex, edge::Hedge)
Add the edge to the list of edges leaving the vert.
"""
function setleavingedge!(vert::Vertex, edge::Hedge)
    for i in 1:3
        if ismissing(vert.leavingEdges[i])
            vert.leavingEdges[i] = edge
            break
        end
    end
end # setleavingedge!
export setleavingedge!

"""
    setleavingedge!(vert::Vertex, oldedge::Hedge, newedge::Hedge)
If two edges are provided, replace old edge by the new in the list of leaving edges.
"""
function setleavingedge!(vert::Vertex, oldedge::Hedge, newedge::Hedge)
    for i in 1:3
        if ismissing(vert.leavingEdges[i]) break end
        if vert.leavingEdges[i] == oldedge
            vert.leavingEdges[i] = newedge
            break
        end
    end
end
export setleavingedge!

function setleavingedges!(cells::Vector{Cell})
    for i in 1:length(cells)
        cell = cells[i]
        for edge in cell
            p1, p2 = edge.originVertex, edge.nextEdge.originVertex
            setleavingedge!(p1, edge)
        end
    end
end
export setleavingedges!

"""
    findthisvert(vert::Vertex, vertices::Vector{Vertex})
Find a vertex in a list of vertices and return its index, and 0 if not found
"""
function findthisvert(vert::Vertex, vertices::Vector{Vertex})
    found = 0 # inde
    for i in 1:length(vertices)
        if vert == vertices[i]
            found = i
            break
        end
    end
    return found
end # findthisvert

"""
    getcellverts_table(lines::Vector{String}, info::Vector{Int})
Builds cell verts table from the processed input file.

The vert reference is the index on the coordinate table.
"""
function getcellverts_table(lines::Vector{String}, info::Vector{Int})
    cellverts = Vector{Vector{Int}}()
    for i in 1:info[1]
        push!(cellverts, [parse(Int,str) for str in split(lines[i])])
    end
    return cellverts
end # getcellverts_table

"""
    getcellverts(cell::Cell)
Transverse a cell and return the list of vertex objects
"""
function getcellverts(cell::Cell)
    verts = Vector{Vertex}()
    for edge in cell
        push!(verts,edge.originVertex)
    end
    return verts
end # getcellverts

"""
    getcellsizes(cellverts::Vector{Vector{Int}})
Gets the number of verts in each cell (cell sizes) from the cell verts table
"""
function getcellsizes(cellverts::Vector{Vector{Int}})
    cellsizes = Vector{Int64}(undef,size(cellverts))
    for i in 1:length(cellverts)
        cellsizes[i] = length(cellverts[i])
    end
    return cellsizes
end # getcellsizes

function getcellsizes(cellverts::Array{Int64,2})
    cellsizes = Vector{Int64}(undef,size(cellverts, 1))
    row_size = size(cellverts, 2)
    for i in 1:size(cellverts, 1)
        cellsize = 0
        for j in 1:row_size
            cellverts[i,j] <= 0 && break
            cellsize += 1
        end
        cellsizes[i] = cellsize
    end
    return cellsizes
end # getcellsizes

"""
    updatecell!(cell::Cell)
Update cell area, perimeter and centroid
"""
function updatecell!(cell::Cell)
    sumArea = 0.0
    sumPerim = 0.0
    xcent = 0.0
    ycent = 0.0

    # loop over vertex and get the measurements
    for edge in cell
        p1, p2 = edge.originVertex, edge.nextEdge.originVertex
        areaterm = p1.x*p2.y - p2.x*p1.y
        sumArea += areaterm
        xcent += (p1.x + p2.x) * areaterm
        ycent += (p1.y + p2.y) * areaterm
        sumPerim += edge.edgeLen
    end

    # Update the values
    cell.areaCell = 0.5*sumArea
    cell.centroid.x = xcent/(6.0*cell.areaCell)
    cell.centroid.y = ycent/(6.0*cell.areaCell)
    cell.perimCell = sumPerim
end # updatecell
export updatecell!

"""
    invertcell!(cell::Cell)
Invert the order of the vertices in a cell object.
"""
function invertcell!(cell::Cell)
    for edge in cell
        edge.nextEdge, edge.prevEdge = edge.prevEdge, edge.nextEdge
        edge = edge.prevEdge
    end
end # invertcell

"""
    getedgeverts(cellverts::Vector{Vector{Int64}})
Get the edge verts table from the cell verts vector of vectors.
"""
function getedgeverts(cellverts::Vector{Vector{Int64}})
    edgeverts = Vector{Vector{Int64}}()
    for i in 1:size(cellverts,1)
        for j in 1:size(cellverts[i],1)-1
            push!(edgeverts,[cellverts[i][j], cellverts[i][j+1]])
        end
        push!(edgeverts,[cellverts[i][size(cellverts[i],1)], cellverts[i][1]])
    end
    return edgeverts
end # getedgeverts

"""
    getedgeverts(cellverts::Array{Int64,2})
Get the edge verts table from the cell verts table (as a matrix).
"""
function getedgeverts(cellverts::Array{Int64,2})
    edgeverts = Vector{Vector{Int64}}()
    for i in 1:size(cellverts,1)
        for j in 1:size(cellverts[i,:],1)-1
            push!(edgeverts,[cellverts[i,j], cellverts[i,j+1]])
        end
        push!(edgeverts,[cellverts[i,length(cellverts[i,:])], cellverts[i,:1]])
    end
    return edgeverts
end # getedgeverts

"""
    findorigin!(edges::Vector{Hedge}, vertices::Vector{Vertex}, vertsInEdge)
Create list of edges and assign the origin vertex to all of them.
"""
function findorigin!(edges::Vector{Hedge}, vertices::Vector{Vertex}, vertsInEdge)
    for i in 1:size(vertsInEdge,1)
        vert = vertices[vertsInEdge[i,1]][1]
        edges[i] = Hedge()
        edges[i].id = i
        edges[i].originVertex = vert
    end
end # findorigin

"""
    findnext!(edges::Vector{Hedge}, vertsinedge, cellsizes)
Assign the next edge field in list of edges.

Requires a list of the verts in all edges and a list pf number of verts per cell.
"""
function findnext!(edges::Vector{Hedge}, vertsInEdge::Vector{Vector{Int64}}, cellsizes::Vector{Int64})
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

"""
    findprevious!(edges::Vector{Hedge}, vertsInEdge, cellsizes)
Assign the previous edge field for each edge in the list provided.

Requires an array with the vert index of each edge and one with the number of vertices per cell.
"""
function findprevious!(edges::Vector{Hedge}, vertsInEdge::Vector{Vector{Int64}}, cellsizes::Vector{Int64})
    edgecell = 1 # count the number of edges that I intereact in a given cell
    icell = length(cellsizes)  # count the cells
    for i in reverse(1:size(vertsInEdge,1))
        if edgecell == cellsizes[icell]
            edges[i].prevEdge = edges[i+cellsizes[icell]-1]
            icell -= 1
            edgecell = 1
            continue
        end
        edges[i].prevEdge = edges[i-1]
        edgecell += 1
    end
end # findnext

"""
    findtwinedges!(edges::Vector{Hedge}, vertsInEdge)
Assign the edge twin field for each edge in the list provided.

Requires an array of vertex index in each edge.
"""
function findtwinedges!(provEdges::Vector{Hedge}, vertsInEdge::Vector{Vector{Int64}})
    ntwins = 0
    for i in 1:length(vertsInEdge)
        verts = reverse(vertsInEdge[i,:][1])
        for j in 1:length(vertsInEdge)
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

"""
    findcontainercell!(edges::Vector{Hedge}, cells::Vector{Cell}, cellsizes)
Assign the container cell field for each edge in the list provided.

Requires an array of vertex index in each edge.
"""
function findcontainercell!(edges::Vector{Hedge}, cells::Vector{Cell}, cellsizes)
    edgeid = 1
    cellid = 1
    for i in 1:size(edges,1)
        if edgeid > cellsizes[cellid]
            edgeid = 1
            cellid += 1
        end
        edges[i].containCell = cells[cellid]
        edgeid += 1
    end
end

"""
    setedgeborder!(edge::Hedge)
Assign the border field for an edge at the border.

Use this run after all the other edge fields had been assigned. A border edge does not have twin!
"""
function setedgeborder!(edge::Hedge)
    edge.border = !isdefined(edge, :twinEdge)
end # assignborders

"""
    newedgelen!(edge::Hedge)
Update edge length
"""
function newedgelen!(edge::Hedge)
    p1 = edge.originVertex
    p2 = edge.nextEdge.originVertex
    edge.edgeLen = distvertices(p1,p2)
end # newedgelen
export newedgelen!

"""
    update_esges!(edge::Hedge, twin::Hedge)
Update pair of edges.
"""
function update_edges!(edge::Hedge, twin::Hedge)
    newedgelen!(edge)
    newedgelen!(twin)
end
export update_edges!

function update_edges!(edge::Hedge, twin::Missing)
    newedgelen!(edge)
end
export update_edges!


"""
    createlistcell(edges::Vector{Hedge}, cellsizes)
Create the list of cells and assign the incident edge for each cell
"""
function createlistcell(edges::Vector{Hedge}, cellsizes)
    cells =  Array{Cell}(undef,size(cellsizes,1))
    edgeid = 1
    for i in 1:size(cellsizes,1)
        cells[i] = Cell()
        cells[i].id = i
        cells[i].incEdge = edges[edgeid]
        cells[i].celltype = 1
        cells[i].centroid = Vertex([0.0,0.0])
        edgeid += cellsizes[i]
    end
    return cells
end # createlistcell

function addtomesh!(mesh::Mesh, vert::Vertex)
    push!(mesh.vertices, vert)
    mesh.lastid[1] += 1
    vert.id = mesh.lastid[1]
end
export addtomesh!

function addtomesh!(mesh::Mesh, edge::Hedge)
    push!(mesh.edges, edge)
    mesh.lastid[2] += 1
    edge.id = mesh.lastid[2]
end
export addtomesh!

function addtomesh!(mesh::Mesh, cell::Cell)
    push!(mesh.cells, cell)
    mesh.lastid[3] += 1
    cell.id = mesh.lastid[3]
end
export addtomesh!

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
    cellverts = getcellverts_table(allLines,cellListInfo)
    cellsizes = getcellsizes(cellverts)
    vertsInEdge = getedgeverts(cellverts)

    # Create list of edges
    edges = Array{Hedge}(undef,size(vertsInEdge,1))
    findorigin!(edges, vertices, vertsInEdge)
    findnext!(edges, vertsInEdge, cellsizes)
    findprevious!(edges, vertsInEdge, cellsizes)
    findtwinedges!(edges,vertsInEdge)

    # Create list of cells
    cells = createlistcell(edges,cellsizes)

    # connect cells and edges lists
    findcontainercell!(edges,cells,cellsizes)

    # Define the largest id for verts, edges and cells
    lastids = [length(vertices), length(edges), length(cells)]
    # get the edges leaving each vert
    setleavingedges!(cells)

    return Mesh(vertices,edges,cells, lastids)
end # importfromfile

function createmesh(cellverts::Array{Int64,2}, coords::Array{Float64,2})

    vertices = getvertexlist(coords)

    # Create list of edges
    vertsInEdge = getedgeverts(cellverts)
    cellsizes = getcellsizes(cellverts)
    edges = Array{Hedge}(undef,size(vertsInEdge,1))
    findorigin!(edges, vertices, vertsInEdge)
    findnext!(edges, vertsInEdge, cellsizes)
    findprevious!(edges, vertsInEdge, cellsizes)
    findtwinedges!(edges,vertsInEdge)

    # Create list of cells
    cells = createlistcell(edges,cellsizes)

    # connect cells and edges lists
    findcontainercell!(edges,cells,cellsizes)

    # Define the largest id for verts, edges and cells
    lastids = [length(vertices), length(edges), length(cells)]
    # get the edges leaving each vert
    setleavingedges!(cells)

    return Mesh(vertices,edges,cells, lastids)
end
export createmesh

# Update the measures of the mesh, namely cell areas and perimeters and edge lens
# Arguments: DCEL object
# Return: DCEL object
function updatesystem!(system::Mesh)

    # Update cells
    for i in 1:length(system.cells)
        updatecell!(system.cells[i])
    end

    # Update edges
    for i in 1:length(system.edges)
        newedgelen!(system.edges[i])
        setedgeborder!(system.edges[i])
    end

    # Update vertices
end # updatesystem

# Export Dcel object to file
# Arguments: file path
# Return: DCEL object
function exporttofile(mesh::Mesh)

    # Create the coordinates list
    coords = extract_vertcoords(mesh)

    # Create the connectivity list
    celltab = vertsincells(mesh)

    open("/home/jhon/Documents/Projects/vertexModelJulia/results/epi.txt", "w") do out
        nvert = length(mesh.vertices)
        write(out, "$nvert\n")
        writedlm(out, coords)
        ncell = length(mesh.cells)
        write(out, "$ncell\n")
        writedlm(out, celltab)
    end
end # exporttofile
export exporttofile

function extract_vertcoords(mesh::Mesh)
    coords = Array{Float64}(undef, length(mesh.vertices), 2)
    for i in 1:length(mesh.vertices)
        coords[i,:] = [mesh.vertices[i].x, mesh.vertices[i].y]
    end
    return coords
end
export extract_vertcoords

function vertsincells(mesh::Mesh)
    celltab = []
    for cell in 1:length(mesh.cells)
        verts_cells = getcellverts(mesh.cells[cell])
        vertids = [findthisvert(verts_cells[vert], mesh.vertices) for vert in 1:length(verts_cells)]
        push!(celltab, vertids)
    end
    return celltab
end
export vertsincells

end # module
