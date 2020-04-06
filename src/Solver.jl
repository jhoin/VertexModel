# This module implements the solver functions
module Solver

# Implements a forward Euler integrator
using ..Meshes
using ..TopologyChanges

# Constants
 const step_size = 0.001
 const dx, dy = 0.0001, 0.0001
 const directions = (90.0*π/360.0, 180.0*π/360.0, 270.0*π/360.0, 0.0*π/360.0)

# Update the local cell
# Arguments: Vertex object and the treatment for boundary (nothing)
# Return: A float representing the vertex energy
function uptade_local!(vert::Vertex, edges::Vector{Hedge}, boundary::Mobile)
    for i in 1:length(edges)
        newedgelen!(edges[i])
        !edges[i].border && newedgelen!(edges[i].twinEdge)
        updatecell!(edges[i].containCell)
    end
end

function uptade_local!(vert::Vertex, edges::Vector{Hedge}, boundary::Spring)
    for i in 1:length(edges)
        newedgelen!(edges[i])
        !edges[i].border && newedgelen!(edges[i].twinEdge)
        updatecell!(edges[i].containCell)
        update_spring!(boundary)
    end
end

function uptade_local!(mesh::SubMesh, boundary::Mobile)
    for i in 1:length(mesh)
        #newedgelen!(mesh.edges[i])
        #!mesh.edges[i].border && newedgelen!(mesh.edges[i].twinEdge)
        update_edges!(mesh.edges[i], mesh.edges[i].twinEdge)
        updatecell!(mesh.cells[i])
    end
end

function uptade_local!(mesh::SubMesh, boundary::Spring)
    for i in 1:length(mesh)
        #newedgelen!(mesh.edges[i])
        #!mesh.edges[i].border && newedgelen!(mesh.edges[i].twinEdge)
        update_edges!(mesh.edges[i], mesh.edges[i].twinEdge)
        updatecell!(mesh.cells[i])
        update_spring!(boundary)
    end
end

# Get a missing energy value for the fixed vertices
# Arguments: vertex object and the treatment boundary
# Return: Missing
function energyvert!(energy::Vector{Float64}, mesh::SubMesh, boundary::Immobile)
    return Missing
end

function energyvert!(energy::Vector{Float64}, edges::Vector{Hedge}, vert::Vertex, boundary::Mobile, K::Float64)
    fill!(energy, 0.0)
    for i in 1:length(edges)
        edge = edges[i]
        cell = edge.containCell
        energy[1] = energy[1] + 0.5*(cell.areaCell/K - 1.0)^2
        energy[2] = energy[2] + 0.5*cell.perim_elast*cell.perimCell^2 / K
        energy[3] = energy[3] + edge.eqbond*edge.edgeLen/sqrt(K)
    end
    total_energy = sum(energy)
    return total_energy
end # energyvert

function energyvert!(energy::Vector{Float64}, mesh::SubMesh, boundary::Spring)
    fill!(energy, 0.0)
    K = mesh.unitlen
    for i in 1:length(mesh)
        edge = mesh.edges[i]
        cell = mesh.cells[i]
        energy[1] = energy[1] + 0.5*(cell.areaCell/K - 1.0)^2
        energy[2] = energy[2] + 0.5*cell.perim_elast*cell.perimCell^2 / K
        energy[3] = energy[3] + edge.eqbond*edge.edgeLen/sqrt(K)
        energy[4] = energy[4] + boundary.spring_const * boundary.len / sqrt(K)
    end
    return sum(energy)
end # energyvert

function energyvert!(energy::Vector{Float64}, mesh::SubMesh, boundary::Mobile)
    fill!(energy, 0.0)
    K = mesh.unitlen
    for i in 1:length(mesh)
        edge = mesh.edges[i]
        cell = mesh.cells[i]
        energy[1] = energy[1] + 0.5*(cell.areaCell/K - 1.0)^2
        energy[2] = energy[2] + 0.5*cell.perim_elast*cell.perimCell^2 / K
        energy[3] = energy[3] + edge.eqbond*edge.edgeLen/sqrt(K)
    end
    return sum(energy)
end # energyvert

function energyvert!(energy::Vector{Float64}, edges::Vector{Hedge}, vert::Vertex, boundary::Spring, K::Float64)
    fill!(energy, 0.0)
    for i in 1:length(edges)
        edge = edges[i]
        cell = edge.containCell
        energy[1] = energy[1] + 0.5(cell.areaCell/K - 1.0)^2
        energy[2] = energy[2] + 0.5*cell.perim_elast*cell.perimCell^2 / K
        energy[3] = energy[3] + edge.eqbond*edge.edgeLen/sqrt(K)
        energy[4] = energy[4] + boundary.spring_const * boundary.len / sqrt(K)
    end
    total_energy = sum(energy)
    return total_energy
end # energyvert

"""
    solve!(mesh::Mesh, t_final::Float64, K::Float64)
Solve the system of equations using the forward Euler method. Updates the mesh type
"""
function solve!(mesh::Mesh, t_final::Float64, K::Float64)
    n_iter = 0
    t = 0.0
    energy_terms = zeros(Float64, (4))
    newmesh = mesh
    while(t < t_final)
        t = n_iter*step_size
        update_topology!(newmesh, 10., K, K*0.1)
        @inbounds for i in 1:length(newmesh.vertices)
            vert = newmesh.vertices[i]
            submesh = get_submesh(vert, K)
            old_energy = energyvert!(energy_terms, submesh, vert.treatBoundary)
            ismissing(old_energy) && continue
            vert_displace!(energy_terms, submesh, old_energy)
        end
        n_iter += 1
        mesh = newmesh
    end
end #solve!
export solve!

"""
    solve!(mesh::Mesh, n_iter::Int64, K::Float64)
Solve the equations using the forward Euler method.
This uses the number of iterations instead of final time.
"""
function solve!(mesh::Mesh, n_iter::Int64, K::Float64)
    iter = 0
    t = 0.0
    energy_terms = zeros(Float64, (4))
    newmesh = mesh
    while(iter < n_iter)
        t = iter*step_size
        update_topology!(newmesh, 10., K, K*0.1)
        for i in 1:length(newmesh.vertices)
            vert = newmesh.vertices[i]
            submesh = get_submesh(vert, K)
            old_energy = energyvert!(energy_terms, submesh, vert.treatBoundary)
            ismissing(old_energy) && continue
            vert_displace!(energy_terms, submesh, old_energy)
        end
        iter += 1
        mesh = newmesh
    end
end #solve!
export solve!

function vert_displace!(energy::Vector{Float64},vert::Vertex, d::Float64, old_energy::Float64, K::Float64)
    directions = [90.0*π/360.0, 180.0*π/360.0, 270.0*π/360.0, 0.0*π/360.0]
    edges = collect(skipmissing(vert.leavingEdges))
    coords = vert.x, vert.y
    for i in 1:4
        vert.x = vert.x + d*cos(directions[i])
        vert.y = vert.y + d*sin(directions[i])
        uptade_local!(vert, edges, vert.treatBoundary)
        new_energy = energyvert!(energy, edges, vert, vert.treatBoundary, K)
        new_energy < old_energy && break
        vert.x, vert.y = coords
    end
end

function vert_displace!(energy::Vector{Float64}, mesh::SubMesh, old_energy::Float64)
    vert = mesh.vert
    coords = vert.x, vert.y
    for i in 1:4
        vert.x = vert.x + dx*cos(directions[i])
        vert.y = vert.y + dx*sin(directions[i])
        uptade_local!(mesh, vert.treatBoundary)
        new_energy = energyvert!(energy, mesh, vert.treatBoundary)
        new_energy > old_energy && continue
        vert.x, vert.y = coords
    end
end

end # module
