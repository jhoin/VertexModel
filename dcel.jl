# This module contains the doubly-connected edge list data structure
# It is a provisory file, each object will be moved to its respective file in the near future

module DCEL

    # This object holds the cell information
    mutable struct Cell
        incEdge::Hedge
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
    mutable struct Vertex
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
    function importfromfile(filePath::String)
        open(filePath) do inFile

            # Get the point coordinates
            nPoints = parse(Int,readline(inFile))
            pointCoords = Array{Float64}(undef,2)
            for i in 1:nPoints
                pointCoords[:] = [parse(Float64,str) for str in split(readline(inFile))]
                vertex = Vertex()
                vertex.x = pointCoords[1]
                vertex.y = pointCoords[2]
                vertices[i] = vertex
            end

            # Get the element connectivity matrix
            cellListInfo = [parse(Int,str) for str in split(readline(inFile))]
            vertsInCells = Matrix{Int}(undef,cellListInfo[1],cellListInfo[2])
            for i in 1:cellListInfo[1]
                vertsInCells[i,:] = [parse(Int,str) for str in split(readline(inFile))]
            end

            # Create the vertex list
            # since it is in an array, it makes easier to search for the vertex idexes later
            vertices = Vector{Vertex}(undef,nPoints)
            for i in 1:cellListInfo[1]

            end
                #Vertex(pointCoords[1],pointCoords[2],undef)

        end
    end

    importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/tests/ex.points")
end    

# TODO:
# Import mesh from file
# Update mesh measurements (area, perimeter, lenght...)
# Export mesh to file
# Plotting operations for mesh
