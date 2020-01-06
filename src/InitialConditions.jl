module InitialConditions

using ..Meshes

"""
    rosetta_hexcenters(cellradius::Float64)
Create center points mesh (make figures)
"""
function rosetta_hexcenters(cellradius::Float64)
    w = sqrt(3.) * cellradius # cell width
    h = 2. * cellradius # cell height
    xy = [0.0, 0.0]
    centers = Array{Float64,2}(undef, 8, 2)
    centers[1,:] = [1.5w + 20., 1.75h + 20.]
    dx = [-0.5w, 0.5w, w, 0.5w, -0.5w, -w]
    dy = [-0.75h, -0.75h, 0.0, 0.75h, 0.75h, 0.0]
    for i in 1:6
        centers[i+1,1] = centers[1,1] + dx[i]
        centers[i+1,2] = centers[1,2] + dy[i]
    end
    return centers
end

function rosetta_mesh(cellradius::Float64)
    centers = rosetta_hexcenters(cellradius)
    ncells = 7
    hex_coords = hextable_coords(centers, ncells, cellradius)
    vert_coords = verttable(centers, ncells, hex_coords)
    hex_verts = hextable(vert_coords, ncells, hex_coords)
    return createmesh(hex_verts, vert_coords)
end
export rosetta_mesh

"""
    hexagonal_mesh(dims::Vector{Int64}, rcell::Float64)
Create a hexagonal tiling, with dims[1]xdims[2] cells of radius rcell
"""
function honeycomb_mesh(dims::Vector{Int64}, cellradius::Float64)
    centers = hexagon_centers(dims, cellradius)
    ncells = dims[1]*dims[2]
    hex_coords = hextable_coords(centers, ncells, cellradius)
    vert_coords = verttable(centers, ncells, hex_coords)
    hex_verts = hextable(vert_coords, ncells, hex_coords)
    return createmesh(hex_verts, vert_coords)
end
export honeycomb_mesh

function hextable_coords(centers::Array{Float64,2}, ncells::Int64, cellradius::Float64)
    hex_coords = Array{Float64,3}(undef, ncells, 6, 2)
    for cell in 1:ncells
        hex = draw_hex(centers[cell,:], cellradius)
        hex_coords[cell, :, :] = hex
    end
    return hex_coords
end

function verttable(centers::Array{Float64,2}, ncells::Int64, hex::Array{Float64,3})
    verts = Array{Float64}(undef, ncells*6, 2)
    cur_vert = 1
    for cell in 1:ncells
        verts[cur_vert:cur_vert + 5, :] = hex[cell, :, :]
        cur_vert += 6
    end
    return unique(verts, dims=1)
end

function hextable(verts::Array{Float64,2}, ncells::Int64, hex_coords::Array{Float64,3})
    hex_mesh = Array{Int64}(undef, ncells, 6)
    for cell in 1:ncells
        for vert in 1:6
            vert_id = find_coord(hex_coords[cell, vert, :], verts)
            hex_mesh[cell, vert] = vert_id
        end
    end
    return hex_mesh
end

function hexagon_centers(dims::Vector{Int64}, cellradius::Float64 )
    w = sqrt(3.) * cellradius # cell width
    h = 2. * cellradius # cell height
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

function draw_hex(c::Vector{Float64}, radius::Float64)
    hex_coords = Array{Float64}(undef, 6, 2)
    directions = [6., 5., 4., 3., 2., 1.]
    for i in 1:length(directions)
        angle_deg = 60.0 * i + 30.0
        ang = π/180.0 * angle_deg
        hex_coords[i, :] = [c[1] + radius*cos(ang), c[2] + radius*sin(ang)]
    end
    return hex_coords
end

function find_coord(coord::Vector{Float64}, coords::Array{Float64,2})
    found = 0
    for vert in 1:size(coords, 1)
        if coord ≈ coords[vert,:]
            found = vert
            break
        end
    end
    return found
end

end #module
