# This module implements the solver functions
module Solver

# Implements a forward Euler integrator
using ..Meshes
using ..TopologyChanges

# Update the local cell
# Arguments: Vertex object and the treatment for boundary (nothing)
# Return: A float representing the vertex energy
function uptade_local!(vert::Vertex)
    for i in 1:3
        if !isassigned(vert.leavingEdges, i) break end
        edge = vert.leavingEdges[i]
        cell = edge.containCell
        newedgelen!(edge)
        !edge.border && newedgelen!(edge.twinEdge)
        updatecell!(cell)
        typeof(vert.treatBoundary) <: Spring && update_spring!(vert)
    end
end

# Calculate the energy of a vertex for vertices without boundary treatment
# Arguments: Vertex object and the treatment for boundary (nothing)
# Return: A float representing the vertex energy
function energyvert(vert::Vertex, boundary::Mobile)
    energy_area = 0.0
    energy_perim = 0.0
    bond_energy = 0.0
    for i in 1:3
        if !isassigned(vert.leavingEdges, i) break end
        edge = vert.leavingEdges[i]
        cell = edge.containCell
        energy_area += 0.5*cell.area_elast*(cell.areaCell - cell.eqarea)^2
        energy_perim += 0.5*cell.perim_elast*cell.perimCell^2
        bond_energy += edge.eqbond*edge.edgeLen
    end
    total_energy = energy_area + energy_perim + bond_energy
    return total_energy
end # energyvert

# Get a missing energy value for the fixed vertices
# Arguments: vertex object and the treatment boundary
# Return: Missing
function energyvert(vert::Vertex, voundary::Immobile)
    return Missing
end

# Calculate the energy of a vertex for vertices with spring boundary condition
# Arguments: Vertex object and the treatment for boundary (spring type)
# Return: A float representing the vertex energy
function energyvert(vert::Vertex, boundary::Spring)
    energy_area = 0.0
    energy_perim = 0.0
    bond_energy = 0.0
    spring_term = 0.
    for i in 1:3
        if !isassigned(vert.leavingEdges, i) break end
        edge = vert.leavingEdges[i]
        cell = edge.containCell
        energy_area += 0.5*cell.area_elast*(cell.areaCell - cell.eqarea)^2
        spring_term = -vert.treatBoundary.spring_const * vert.treatBoundary.displaced
        energy_perim += 0.5*cell.perim_elast*cell.perimCell^2
        bond_energy += edge.eqbond*edge.edgeLen
    end
    total_energy = energy_area + energy_perim + bond_energy + spring_term
    return total_energy
end # energyvert

function energyvert(vert::Vertex, boundary::Spring, K::Float64)
    energy_area = 0.0
    energy_perim = 0.0
    bond_energy = 0.0
    spring_term = 0.
    rootK = sqrt(K)
    for i in 1:3
        if !isassigned(vert.leavingEdges, i) break end
        edge = vert.leavingEdges[i]
        cell = edge.containCell
        energy_area += 0.5*(cell.areaCell/K - 1.0)^2
        spring_term = vert.treatBoundary.spring_const*vert.treatBoundary.len / rootK
        energy_perim += 0.5*cell.perim_elast*cell.perimCell^2 / K
        bond_energy += edge.eqbond*edge.edgeLen/rootK
    end
    total_energy = energy_area + energy_perim + bond_energy + spring_term
    return total_energy
end # energyvert

function energyvert(vert::Vertex, boundary::Mobile, K::Float64)
    energy_area = 0.0
    energy_perim = 0.0
    bond_energy = 0.0
    spring_term = 0.
    rootK = sqrt(K)
    for i in 1:3
        if !isassigned(vert.leavingEdges, i) break end
        edge = vert.leavingEdges[i]
        cell = edge.containCell
        energy_area += 0.5*(cell.areaCell/K - 1.0)^2
        energy_perim += 0.5*cell.perim_elast*cell.perimCell^2 / K
        bond_energy += edge.eqbond*edge.edgeLen/rootK
    end
    total_energy = energy_area + energy_perim + bond_energy + spring_term
    return total_energy
end # energyvert

# Solve the system of equations using the forward Euler method
# Arguments: Mesh object and final time
# Return: Mesh object updated
function solve!(mesh::Mesh, t_final::Float64)
    n_iter = 0
    step_size = 0.001
    t = 0.0
    newmesh = mesh
    while(t < t_final)
        t = n_iter*step_size
        #update_topology!(newmesh, 2.)

        # Update vertex positions
        dx, dy = 0.001, 0.001
        for i in 1:length(newmesh.vertices)
            vert = newmesh.vertices[i]
            K = vert.leavingEdges[1].containCell.eqarea
            old_energy = energyvert(mesh.vertices[i], mesh.vertices[i].treatBoundary, K)
            ismissing(old_energy) && continue
            vert_displace!(vert, dx, old_energy, K)
        end
        n_iter += 1
        mesh = newmesh
    end
    return mesh
end #solve!
export solve!

function solve!(mesh::Mesh, n_iter::Int64)
    iter = 1
    while iter <= n_iter

        # Update all new variables!
        newmesh = mesh
        update_topology!(newmesh, 2.0)

        # Get a random vertex and calculate its energy
        vert = newmesh.vertices[1 + trunc(Int, (length(newmesh.vertices)+1-1)*rand())]
        K = vert.leavingEdges[1].containCell.eqarea
        old_energy = energyvert(vert, vert.treatBoundary, K)

        # Move the vertex
        attempt_move!(vert, 0.3)
        uptade_local!(vert)

        # Calculate the vertex energy after moving
        new_energy = energyvert(vert, vert.treatBoundary, K)

        # Accept this iteration?
        #isAccepted = "N"
        #call random_number(rnd)
        #probAccept = exp((newVertEnergy - oldVertEnergy)/T)
        #call random_number(rnd)
        probAccept = 0.05
        if new_energy < old_energy || probAccept < rand()

            #Update all variables!
            mesh = newmesh

            #Mark the iteration as accepted
            #isAccepted = "Y"
        end
        iter = iter + 1
    end
end
export solve!

function attempt_move!(vert::Vertex, radius::Float64)
    rnd1 = rand()
    rnd2 = rand()

    # If one is bigger than the other, swap them
    if rnd2 > rnd1 rnd2, rnd1 = rnd1, rnd2 end

    # Moooooove the points
    angle = 360*rand()
    vert.x = vert.x + rnd2*radius*cos(angle*π*180.0)
    vert.y = vert.y + rnd2*radius*sin(angle*π*180.0)
    #vert.x = vert.x + rnd2*radius*cos(2.0*π*rnd1/rnd2)
    #vert.y = vert.y + rnd2*radius*sin(2.0*π*rnd1/rnd2)
end

function vert_displace!(vert::Vertex, d::Float64, old_energy::Float64, K::Float64)
    directions = [90.0*π/360.0, 180.0*π/360.0, 270.0*π/360.0, 0.0*π/360.0]
    coords = vert.x, vert.y

    for i in 1:4
        vert.x = vert.x + d*cos(directions[i])
        vert.y = vert.y + d*sin(directions[i])
        uptade_local!(vert)
        new_energy = energyvert(vert, vert.treatBoundary, K)
        new_energy < old_energy && break
        vert.x, vert.y = coords
    end
end

# Calculate the vertex displacement along the x axis
# Arguments: vertex, an infinitesimal displacement and the vert energy before displacement
# Return: displacement along the x axis
function vert_displacex(vert::Vertex, dx::Float64, old_energy::Float64, K)
    x = vert.x
    vert.x = vert.x + dx
    uptade_local!(vert)
    dx_plus = energyvert(vert, vert.treatBoundary, K)
    vert.x = x

    x = vert.x
    vert.x = vert.x - dx
    uptade_local!(vert)
    dx_minus = energyvert(vert, vert.treatBoundary, K)
    vert.x = vert.x + dx
    vert.x = x

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
function vert_displacey(vert::Vertex, dy::Float64, old_energy::Float64, K)
    y = vert.y
    vert.y = vert.y + dy
    uptade_local!(vert)
    dy_plus = energyvert(vert, vert.treatBoundary, K)
    vert.y = vert.y - dy
    vert.y = y

    y = vert.y
    vert.y = vert.y - dy
    uptade_local!(vert)
    dy_minus = energyvert(vert, vert.treatBoundary, K)
    vert.y = vert.y + dy
    vert.y = y

    # Update or not the positions
    move_sign = one(1.0)
    if(dy_plus > dy_minus)
        move_sign = -one(1.0)
    end
    return move_sign
end # vert_displacey

end # module
