
module Solver

# Implements a forward Euler integrator
using ..DCEL
using Traceur
function uptade_local!(vert::Vertex)
    for i in 1:length(vert.leavingEdges)
        if !isassigned(vert.leavingEdges, i) break end
        edge = vert.leavingEdges[i]
        cell = edge.containCell
        newedgelen!(edge)
        updatecell!(cell)
    end
end

# Calculate the energy of a vertex
# Arguments: Vertex object
# Return: A float representing the vertex energy
function energyvert(vert::Vertex)
    energy_area = 0.0
    energy_perim = 0.0
    bond_energy = 0.0
    for i in 1:length(vert.leavingEdges)
        if !isassigned(vert.leavingEdges, i) break end
        edge = vert.leavingEdges[i]
        cell = edge.containCell
        energy_area += 0.5*cell.area_elast*(cell.areaCell - cell.eqarea)^2
        energy_perim += 0.5*cell.perim_elast*(cell.perimCell - cell.eqperim)^2
        bond_energy += edge.eqbond*edge.edgeLen
    end
    total_energy = energy_area + energy_perim + bond_energy
    return total_energy
end # energyvert

# Solve the system of equations using the forward Euler method
# Arguments: Mesh object and final time
# Return: Mesh object updated
function solve!(mesh::Dcel, t_final::Float64)
    n_iter = 0
    step_size = 0.001
    t = 0.0
    newmesh = mesh
    while(t < t_final)
        t = n_iter*step_size

        # Update vertex positions
        dx, dy = 0.0001, 0.0001
        for i in 1:length(newmesh.listVert)
            vert = newmesh.listVert[i]
            old_energy = energyvert(mesh.listVert[i])
            disx = vert_displacex(vert, dx, old_energy)
            disy = vert_displacey(vert, dy, old_energy)
            vert.x = vert.x + disx*dx
            vert.y = vert.y + disy*dy
            uptade_local!(vert)
        end
        n_iter += 1
        mesh = newmesh
    end
    return mesh
end #solve!
export solve!

# Calculate the vertex displacement along the x axis
# Arguments: vertex, an infinitesimal displacement and the vert energy before displacement
# Return: displacement along the x axis
function vert_displacex(vert::Vertex, dx::Float64, old_energy::Float64)
    vert.x = vert.x + dx
    uptade_local!(vert)
    dx_plus = energyvert(vert)
    vert.x = vert.x - dx

    vert.x = vert.x - dx
    uptade_local!(vert)
    dx_minus = energyvert(vert)
    vert.x = vert.x + dx

    # Update or not the positions
    move_sign = one(1.0)
    if(dx_plus > dx_minus)
        move_sign = -one(1.0)
    end
    return move_sign
end # vert_displacex

# Calculate the vertex displacement along the y axis
# Arguments: vertex, an infinitesimal displacement and the vert energy before displacement
# Return: displacement along the y axis
function vert_displacey(vert::Vertex, dy::Float64, old_energy::Float64)
    vert.y = vert.y + dy
    uptade_local!(vert)
    dy_plus = energyvert(vert)
    vert.y = vert.y - dy

    vert.y = vert.y - dy
    uptade_local!(vert)
    dy_minus = energyvert(vert)
    vert.y = vert.y + dy

    # Update or not the positions
    move_sign = one(1.0)
    if(dy_plus > dy_minus)
        move_sign = -one(1.0)
    end
    return move_sign
end # vert_displacey

end # module
