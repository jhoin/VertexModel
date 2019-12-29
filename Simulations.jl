module Simulations

using ..Meshes
using ..PointGeometry
using ..InitialConditions
using ..Solver
using Gadfly, Compose

# Initialise model parameters
# Arguments: Mesh object
# Return: Updated mesh
function init_model!(mesh::Mesh)
    for i in 1:length(mesh.cells)
        cell = mesh.cells[i]
        cell.area_elast = 0.5
        cell.perim_elast = 0.05
    end

    for i in 1:length(mesh.edges)
        edge = mesh.edges[i]
        edge.eqbond = 25.0
    end

    boundingbox!(mesh)
end
export init_model!

function apoptosis()
    mesh = rosetta_mesh(30.)
    updatesystem!(mesh)
    adhesion_matrix = [0.2 0.8; 0.8 0.2]
    plot_mesh(mesh)

    # Define cell types
    a_cell = mesh.cells[1]
    a_cell.celltype = 2


    # Paremeters
    A_0 = 3000.
    K = 0.5
    for i in 1:length(mesh.cells)
        mesh.cells[i].perim_elast = 0.5/K*A_0
    end

    for i in 1:length(mesh.edges)
        edge = mesh.edges[i]
        cell = edge.containCell
        if edge.border
            edge.eqbond = 0.5/ K*sqrt(A_0^3)
            continue
        end
        twincell = edge.twinEdge.containCell
        bondenergy = adhesion_matrix[cell.celltype, twincell.celltype]
        edge.eqbond = bondenergy/ K*sqrt(A_0^3)
    end
    boundingbox!(mesh)

    a_vert = a_cell.incEdge.originVertex
    energies = zeros(Float64, (4))
    edges = collect(skipmissing(a_vert.leavingEdges))
    submesh = get_submesh(a_vert, A_0)
    solve!(mesh, 30.0, A_0)
    #println(a_cell.incEdge.edgeLen)
    return mesh
end
export apoptosis

function boundingbox!(mesh::Mesh)

    # The vertices of the box
    coords = extract_vertcoords(mesh)
    xmax = maximum(coords[:,1])*1.1
    xmin = minimum(coords[:,1])*1.1
    ymax = maximum(coords[:,2])*1.1
    ymin = minimum(coords[:,2])*1.1

    # Set the anchring positions and create the spring object
    springs = Vector{Spring}()
    for i in 1:length(mesh.vertices)
        vert = mesh.vertices[i]
        if ismissing(vert.leavingEdges[3])
            anchors = anchortobox(vert, [xmax, xmin, ymax, ymin])
            springlen = sqrt((vert.x - anchors[1])^2 + (vert.y-anchors[2])^2)
            spring = Spring(anchors[1], anchors[2], 0.35, springlen, vert)
            push!(springs, spring)
            vert.treatBoundary = spring
        else
            vert.treatBoundary = Mobile(vert)
        end
    end

    # NOTE: Necesary post processing (MANUAL)
    #vert = springs[9].vert
    #vert.treatBoundary = Mobile(vert)
end

function anchortobox(vert::Vertex, boxxy::Vector{Float64})
    topline = [boxxy[1] boxxy[4]; boxxy[2] boxxy[4]]
    leftline = [boxxy[2] boxxy[3]; boxxy[2] boxxy[4]]
    bottomline = [boxxy[2] boxxy[3]; boxxy[1] boxxy[3]]
    rightline = [boxxy[1] boxxy[4]; boxxy[1] boxxy[3]]

    boxlines = [topline, leftline, bottomline, rightline]
    directions = [π/2., π, 3. *π/2., 2. *π]
    dist = Inf
    coords = [0. 0.]
    for i in 1:4
        x = vert.x + 0.05*vert.x*cos(directions[i])
        y = vert.y + 0.05*vert.y*sin(directions[i])
        coords_new = intersection([vert.x vert.y; x y], boxlines[i])
        dist_new = sqrt((vert.x-coords_new[1])^2 + (vert.y-coords_new[2])^2)
        if dist_new < dist
            coords = coords_new
            dist = dist_new
        end
    end
    return coords
end

# Create a point to anchor the spring
# Arguments: vertex object
# Return: a pair of floats with the new
function anchor_coords(vert::Vertex)
    a_cell = vert.leavingEdges[1].containCell
    an_edge = vert.leavingEdges[1].nextEdge.originVertex
    coords = rotatepoint(vert, an_edge)
    return coords
end

# Plot mesh
# Arguments: mesh
# Return: Plot
function plot_mesh(mesh::Mesh)
    coords = extract_vertcoords(mesh)
    cell_verts = vertsincells(mesh)

    # Fill up x and y coordinates for each cell
    xs = []
    ys = []
    for cell in 1:length(cell_verts)
        ids = cell_verts[cell]
        xs_cell = [coords[ids[vert], 1] for vert in 1:length(ids)]
        ys_cell = [coords[ids[vert], 2] for vert in 1:length(ids)]
        push!(xs, xs_cell)
        push!(ys, ys_cell)
    end

    # plot the grid
    #plt_mesh = plot(x=xs[1], y=ys[1], Geom.polygon(preserve_order=true), Guide.annotation(compose(context(), text(mesh.cells[1].centroid.x, mesh.cells[1].centroid.y,   "$(mesh.cells[1].id)"))))
    plt_mesh = plot(x=xs[1], y=ys[1], Geom.polygon(preserve_order=true))
    for cell in 2:length(xs)
        push!(plt_mesh, layer(x=xs[cell], y=ys[cell], Geom.polygon(preserve_order=true)))
        #push!(plt_mesh, Guide.annotation(compose(context(), text(mesh.cells[cell].centroid.x, mesh.cells[cell].centroid.y, "$(mesh.cells[cell].id)"))))
    end
    #display(plt_mesh)
    #plot_anchors!(mesh, plt_mesh)
    display(plt_mesh)
end
export plot_mesh

function plot_grid(grid)

    plt_mesh = plot(x=grid[1,:,1], y=grid[1,:,2], Geom.polygon(preserve_order=true))
    for cell in 2:length(grid[:,1,1])
        push!(plt_mesh, layer(x=grid[cell,:,1], y=grid[cell,:,2], Geom.polygon(preserve_order=true)))
    end
    display(plt_mesh)
end

function pts(coords::Array{Float64})
    plt_mesh = plot(x=coords[:,1], y=coords[:,2])
    display(plt_mesh)
end

function plot_anchors!(mesh::Mesh, plt)
    xanchors = Vector{Float64}()
    yanchors = Vector{Float64}()
    for i in 1:length(mesh.vertices)
        vert = mesh.vertices[i]
        if isa(vert.treatBoundary, Spring)
            spring = vert.treatBoundary
            push!(xanchors, vert.treatBoundary.x)
            push!(yanchors, vert.treatBoundary.y)
        end
    end
    push!(plt, layer(x=xanchors, y=yanchors))
end

end
