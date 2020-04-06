# main.jl
include("Meshes.jl")
include("InitialConditions.jl")
include("PointGeometry.jl")
include("TopologyChanges.jl")
include("Solver.jl")
include("Simulations.jl")

using .Meshes
using .InitialConditions
using .TopologyChanges
using .Solver
using .Simulations

function main()
    mesh = apoptosis()
end

@time main()

# TODO:
# Fix the leaving edge assignment. None of the topological operations are taking this parameter into account
# Write tests for the T1 topological transistion
# Implement T2 transistion
# Implement periodic boundary conditions
