using Test
using .Meshes
using .TopologyChanges
using .PointGeometry

function test_cell(cell::Cell)
    @testset "Test cell" begin
        for edge in cell
            @test edge.containCell == cell
            @test edge.id == edge.prevEdge.nextEdge.id
            @test edge.prevEdge.containCell.id == cell.id
            @test edge.nextEdge.containCell.id == cell.id
            orig = edge.originVertex
            @test orig != edge.nextEdge.originVertex
            @test orig != edge.prevEdge.originVertex
            @test isleavingedge(orig, edge)
        end
    end
end

@testset "Test all cells" begin
    mesh = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/results/ic.txt")
    updatesystem!(mesh)
    for i in 1:length(mesh.cells)
        cell = mesh.cells[i]
        for edge in cell
            @test edge.containCell == cell
            @test edge.prevEdge.containCell.id == cell.id
            orig = edge.originVertex
            @test isleavingedge(orig, edge)
        end
    end
end

function isleavingedge(vert::Vertex, edge::Hedge)
    isleavingedge = false
    for i in 1:3
        if !isassigned(vert.leavingEdges, i) break end
        if vert.leavingEdges[i] == edge isleavingedge = true end
    end
    return isleavingedge
end

@testset "Test cells after division" begin
    TopologyChanges.cell_division!(mesh.cells[13])
    for i in 1:length(mesh.cells)
        cell = mesh.cells[i]
        for edge in cell
            @test edge.containCell == cell
            @test edge.prevEdge.containCell.id == cell.id
            orig = edge.originVertex
            @test isleavingedge(orig, edge)
        end
    end
end

@testset "test after transition" begin
    mesh = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/results/ic.txt")
    Meshes.updatesystem!(mesh)
    focal_edge = mesh.edges[114]

    TopologyChanges.t1transition!(focal_edge)
    TopologyChanges.t1transition!(foca.twinEdge)
end

@testset "T1 transition" begin
    mesh = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/results/ic.txt")
    Meshes.updatesystem!(mesh)

    focal_edge = mesh.edges[114]
    @test focal_edge != focal_edge.containCell.incEdge
    @test !focal_edge.border
    @test !focal_edge.nextEdge.border
    @test !focal_edge.prevEdge.border
    @test !TopologyChanges.check_t1(focal_edge)
    p1 = focal_edge.originVertex
    p2 = focal_edge.nextEdge.originVertex
    PointGeometry.rotatepoints!(p1, p2)

    # Assert the incident edge field in cell object
    focal_cell = focal_edge.containCell
    if focal_cell.incEdge == focal_edge
        focal_cell.incEdge = focal_edge.nextEdge
    end
    @test focal_edge.containCell == focal_cell

    # Store some info to test against later
    focal_next = focal_edge.nextEdge
    focal_prev = focal_edge.prevEdge

    #@test focal_edge.prevEdge.border
    edge = focal_edge.prevEdge.twinEdge
    adj_cell1 = edge.containCell
    TopologyChanges.edgeremove!(focal_edge, edge)
    @test focal_edge.containCell == edge.containCell
    @test focal_edge.containCell != focal_edge.prevEdge.containCell
    TopologyChanges.addedgeat!(focal_edge, edge)
    @test focal_edge.nextEdge == edge
    @test focal_edge.nextEdge.originVertex == p2
    @test focal_edge.twinEdge.originVertex == focal_edge.nextEdge.originVertex
    @test focal_prev != focal_edge.prevEdge
    @test focal_next != focal_edge.nextEdge

    focal_edge = focal_edge.twinEdge
    @test focal_edge != focal_edge.containCell.incEdge
    edge = focal_edge.prevEdge.twinEdge
    adj_cell2 = edge.containCell
    @test adj_cell1 != adj_cell2
    TopologyChanges.edgeremove!(focal_edge, edge)
    @test focal_edge.containCell == edge.containCell
    @test focal_edge.containCell != focal_edge.prevEdge.containCell
    TopologyChanges.addedgeat!(focal_edge, edge)
    @test focal_edge.nextEdge == edge
    @test focal_edge.nextEdge.originVertex == p1
    @test focal_edge.twinEdge.originVertex == focal_edge.nextEdge.originVertex

    Meshes.setleavingedge!(p1, focal_edge.twinEdge.nextEdge, focal_edge.nextEdge)
    @test isleavingedge( p1, focal_edge.nextEdge)
    focal_edge = focal_edge.twinEdge
    Meshes.setleavingedge!(p2, focal_edge.twinEdge.nextEdge, focal_edge.nextEdge)
    @test isleavingedge( p2, focal_edge.nextEdge)

    # Test cells after transition
    test_cell(focal_edge.containCell)
    test_cell(focal_edge.twinEdge.containCell)
    test_cell(adj_cell1)
    test_cell(adj_cell2)

    ###################### AGAIN #####################

    focal_edge = mesh.edges[114]
    @test focal_edge != focal_edge.containCell.incEdge
    @test !focal_edge.border
    @test !focal_edge.nextEdge.border
    @test !focal_edge.prevEdge.border
    @test !TopologyChanges.check_t1(focal_edge)
    p1 = focal_edge.originVertex
    p2 = focal_edge.nextEdge.originVertex
    PointGeometry.rotatepoints!(p1, p2)

    # Assert the incident edge field in cell object
    focal_cell = focal_edge.containCell
    if focal_cell.incEdge == focal_edge
        focal_cell.incEdge = focal_edge.nextEdge
    end
    @test focal_edge.containCell == focal_cell

    # Store some info to test against later
    focal_next = focal_edge.nextEdge
    focal_prev = focal_edge.prevEdge

    #@test focal_edge.prevEdge.border
    edge = focal_edge.prevEdge.twinEdge
    adj_cell1 = edge.containCell
    TopologyChanges.edgeremove!(focal_edge, edge)
    @test focal_edge.containCell == edge.containCell
    @test focal_edge.containCell != focal_edge.prevEdge.containCell
    TopologyChanges.addedgeat!(focal_edge, edge)
    @test focal_edge.nextEdge == edge
    @test focal_edge.nextEdge.originVertex == p2
    @test focal_edge.twinEdge.originVertex == focal_edge.nextEdge.originVertex
    @test focal_prev != focal_edge.prevEdge
    @test focal_next != focal_edge.nextEdge

    focal_edge = focal_edge.twinEdge
    @test focal_edge != focal_edge.containCell.incEdge
    edge = focal_edge.prevEdge.twinEdge
    adj_cell2 = edge.containCell
    @test adj_cell1 != adj_cell2
    TopologyChanges.edgeremove!(focal_edge, edge)
    @test focal_edge.containCell == edge.containCell
    @test focal_edge.containCell != focal_edge.prevEdge.containCell
    TopologyChanges.addedgeat!(focal_edge, edge)
    @test focal_edge.nextEdge == edge
    @test focal_edge.nextEdge.originVertex == p1
    @test focal_edge.twinEdge.originVertex == focal_edge.nextEdge.originVertex

    Meshes.setleavingedge!(p1, focal_edge.twinEdge.nextEdge, focal_edge.nextEdge)
    @test isleavingedge( p1, focal_edge.nextEdge)
    focal_edge = focal_edge.twinEdge

    Meshes.setleavingedge!(p2, focal_edge.twinEdge.nextEdge, focal_edge.nextEdge)
    @test isleavingedge( p2, focal_edge.nextEdge)

    # Test cells after transition
    test_cell(focal_edge.containCell)
    test_cell(focal_edge.twinEdge.containCell)
    test_cell(adj_cell1)
    test_cell(adj_cell2)


end

@testset "Check T1 transition" begin
    mesh = Meshes.importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/results/ic.txt")
    Meshes.updatesystem!(mesh)

    TopologyChanges.t1transition!(mesh.edges[114], 1.5)

end

@testset "Cell division" begin
    mesh = importfromfile("/home/jhon/Documents/Projects/vertexModelJulia/results/ic.txt")
    updatesystem!(mesh)

    # Test the intersection, by dividing one cell
    cell = mesh.cells[13]
    c2 = shortest_axis(cell)
    @test typeof(c2) <: Vertex
    edges_crossed, verts_crossing = TopologyChanges.intersectcell!(mesh, cell, c2)
    @test length(edges_crossed) == 2
    @test length(verts_crossing) == 2
    @test edges_crossed[1].containCell == cell
    @test edges_crossed[2].containCell == cell

    # Test the split formed by the intersection in the dividing cell
    split1 = TopologyChanges.addedgeat!(mesh, edges_crossed[1], verts_crossing[1])
    @test split1.prevEdge == edges_crossed[1]
    @test edges_crossed[1].nextEdge == split1
    @test split1.containCell == cell
    @test split1.originVertex == verts_crossing[1]
    split2 = TopologyChanges.addedgeat!(mesh, edges_crossed[2], verts_crossing[2])
    @test split2.prevEdge == edges_crossed[2]
    @test edges_crossed[2].nextEdge == split2
    @test split1.containCell == cell
    @test split2.originVertex == verts_crossing[2]

    vert = edges_crossed[2].twinEdge.originVertex

    # Operations to split the twin
    twin_split = TopologyChanges.addedgeat!(mesh, edges_crossed[2].twinEdge.prevEdge)
    twin_split.twinEdge, split2.twinEdge = split2, twin_split
    twin_split.originVertex = edges_crossed[2].twinEdge.originVertex
    edges_crossed[2].twinEdge.originVertex = verts_crossing[2]
    @test twin_split == split2.twinEdge
    @test split2.twinEdge == twin_split
    @test twin_split.containCell.id == twin_split.prevEdge.containCell.id
    @test twin_split.originVertex.x ≈ vert.x

    setleavingedge!(edges_crossed[2].twinEdge.originVertex, edges_crossed[2].twinEdge, twin_split)
    setleavingedge!(vert, edges_crossed[2].twinEdge)
    @test isleavingedge(vert, edges_crossed[2].twinEdge)
end
