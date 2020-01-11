module InitialConditions

using ..Meshes
using ..PointGeometry

"""
    rosetta_hexcenters(rcell::Float64)
Return the center points of each hexagon in a rosetta (make figures)
"""
function rosetta_hexcenters(rcell::Float64)
    w = sqrt(3.) * rcell # cell width
    h = 2. * rcell # cell height
    xy = [0.0, 0.0]
    centers = Array{Float64,2}(undef, 7, 2)
    centers[1,:] = [1.5w + 20., 1.75h + 20.]
    dx = [-0.5w, 0.5w, w, 0.5w, -0.5w, -w]
    dy = [-0.75h, -0.75h, 0.0, 0.75h, 0.75h, 0.0]
    for i in 1:6
        centers[i+1,1] = centers[1,1] + dx[i]
        centers[i+1,2] = centers[1,2] + dy[i]
    end
    return centers
end

"""
    rosetta_mesh(rcell::Float64)
Create a hexagonal rosetta tiling, with one cell at the center surrounded by six other cells.
Each cell has a radius of rcell.
"""
function rosetta_mesh(rcell::Float64)
    centers = rosetta_hexcenters(rcell)
    hex_coords = hextable_coords(centers, rcell)
    vert_coords = verttable(centers, hex_coords)
    hex_verts = hextable(vert_coords, hex_coords)
    return createmesh(hex_verts, vert_coords)
end
export rosetta_mesh

"""
    hexagonal_mesh(dims::Vector{Int64}, rcell::Float64)
Create a hexagonal tiling, with dims[1]xdims[2] cells of radius rcell
"""
function honeycomb_mesh(dims::Vector{Int64}, rcell::Float64)
    centers = hexagon_centers(dims, rcell)
    hex_coords = hextable_coords(centers, rcell)
    vert_coords = verttable(centers, hex_coords)
    hex_verts = hextable(vert_coords, hex_coords)
    return createmesh(hex_verts, vert_coords)
end
export honeycomb_mesh

"""
    hextable_coords(centers::Array{Float64,2}, rcell::Float64)
Create a hexagonal tiling, with dims[1]xdims[2] cells of radius rcell
"""
function hextable_coords(centers::Array{Float64,2}, rcell::Float64)
    n = size(centers,1) # number of hexagons
    hex_coords = Array{Float64,3}(undef, n, 6, 2)
    for cell in 1:n
        hex = draw_hex(centers[cell,:], rcell)
        hex_coords[cell, :, :] = hex
    end
    return hex_coords
end

"""
    verttable(centers::Array{Float64,2}, hex::Array{Float64,3})
Returns a 2D array listing the unique vertices of mesh given the center coordinates of hexagons (centers), and a table storing the coordinates of all vertices of each hexagon.
"""
function verttable(centers::Array{Float64,2}, hex::Array{Float64,3})
    n = size(centers,1) # number of hexagons
    verts = Array{Float64}(undef, 6n, 2)
    cur_vert = 1
    for cell in 1:n
        verts[cur_vert:cur_vert + 5, :] = hex[cell, :, :]
        cur_vert += 6
    end
    return unique(verts, dims=1)
end

"""
    hextable(verts::Array{Float64,2}, hex_coords::Array{Float64,3})
Create a 2D array storing the index of the vertices (from teh matrix verts) presents in each hexagon.
"""
function hextable(verts::Array{Float64,2}, hex_coords::Array{Float64,3})
    n = size(hex_coords, 1) # number of hexagons
    hex_mesh = Array{Int64}(undef, n, 6)
    for cell in 1:n
        for vert in 1:6
            vert_id = find_coord(hex_coords[cell, vert, :], verts)
            hex_mesh[cell, vert] = vert_id
        end
    end
    return hex_mesh
end

"""
    hexagon_centers(dims::Vector{Int64}, rcell::Float64)
Create a 2D array storing the centers of each hexagon.
"""
function hexagon_centers(dims::Vector{Int64}, rcell::Float64)
    w = sqrt(3.) * rcell # cell width
    h = 2. * rcell # cell height
    xy = [0.0, 0.0]
    centers = Array{Float64}(undef, dims[1]*dims[2], 2)
    cell = 1
    for row in 1:dims[1]
        xy[2] += 0.75*h
        for col in 1:dims[2]
            centers[cell,:] = xy .+ 15.
            xy[1] += w
            cell += 1
        end
        xy[1] = 0.5 * w
        if row % 2 == 0
            xy[1] = 0.0
        end
    end
    return centers
end

"""
    draw_hex(c::Vector{Float64}, rcell::Float64)
Create a single hexagon, given the center and radius of the hexagon.
"""
function draw_hex(c::Vector{Float64}, rcell::Float64)
    hex_coords = Array{Float64}(undef, 6, 2)
    directions = [6., 5., 4., 3., 2., 1.]
    for i in 1:length(directions)
        angle_deg = 60.0 * i + 30.0
        ang = π/180.0 * angle_deg
        hex_coords[i, :] = [c[1] + rcell*cos(ang), c[2] + rcell*sin(ang)]
    end
    return hex_coords
end

"""
    find_coord(point::Vector{Float64}, coords::Array{Float64,2})
Find the index position of a point in the list of coordinates. It returns the index if it is found and zero if the point is not there.
"""
function find_coord(point::Vector{Float64}, coords::Array{Float64,2})
    found = 0
    for vert in 1:size(coords, 1)
        if point ≈ coords[vert,:]
            found = vert
            break
        end
    end
    return found
end

"""
    disturb_verts(mesh::Mesh)
Change the coordinates of the vertices in a mesh randomly. Assumes a regular mesh
"""
function disturb_verts!(mesh::Mesh)
    r_disturb = 0.1*mesh.edges[1].edgeLen
    for vert in 1:length(mesh.vertices)
        attempt_move!(mesh.vertices[vert], r_disturb)
    end
end
export disturb_verts!

end #module
