# Hold functions used to test the code!

include("../dcel.jl")
using .DCEL

module TestDCEL

function leavingedge_origvert(mesh)
    # Test leaving edges against origin vertex
    for i in 1:size(mesh.listVert,1)
        vert = mesh.listVert[i]
        for j in 1:3
            if !isassigned(vert.leavingEdges[j]) break end
            edge = vert.leavingEdges[j]
            @assert edge.originVertex == vert "vert leaving edge is different than edge origin vert"
        end
    end
end # leavingedge_origvert
export leavingedge_origvert



end # module

#include("../dcel.jl")
using .TestDCEL

# Test example
mesh = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/results/ex3.points")
updatesystem!(mesh)
leavingedge_origvert(mesh)
