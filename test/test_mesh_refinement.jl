

using Ferrite, Tensors

grid = generate_grid(Triangle, (10,10))
dh = DofHandler(grid)
add!(dh, :u, 3)
close!(dh)

ArrConsNothing = Union{Vector{AffineConstraint{Float64}}, Nothing}


chs = Vector{ArrConsNothing}(nothing, getnnodes(grid))
fildDim = Ferrite.getfielddim(dh, :u)

# chs[1] = [AffineConstraint(1, )]

# check if the node exists

# nodes = grid.nodes
# cells = grid.cells

# refineCell = Vector{Int}
refineCells = [1, 2]

# Triangle element refinement

# newNodes = Vector{Node}[]

function refine_grid_Triangle!(grid::Grid, refineCells::Vector{Int64})
    isirregular = Pair{Int64, Bool}[]
    for cellId in refineCells
        cellCoord = get_cell_coordinates(grid, cellId)
        p1 = (cellCoord[1]+cellCoord[2])/2
        p2 = (cellCoord[1]+cellCoord[3])/2
        p3 = (cellCoord[2]+cellCoord[3])/2
        node1 = Node((p1[1], p1[2]))
        node2 = Node((p2[1], p2[2]))
        node3 = Node((p3[1], p3[2]))
        nnodes = getnnodes(grid)
        id_node1 = nnodes + 1
        id_node2 = nnodes + 2
        id_node3 = nnodes + 3
        push!(grid.cells, Triangle((id_node1, id_node2, id_node3)))
        node1_temp = findall(x->x== node1, nodes)
        if isempty(node1_temp)
            push!(grid.nodes, node1)
            push!(isirregular, id_node1=>true)
        else
            # push!(isirregular, id_node1=>false)
        end
        if isempty(findall(x->x== node2, nodes))
            push!(grid.nodes, node2)
            push!(isirregular, id_node2=>true)
        else
            push!(isirregular, id_node2=>false)
        end

        if isempty(findall(x->x== node3, nodes))
            push!(grid.nodes, node3)
            push!(isirregular, id_node3=>true)
        else
            push!(isirregular, id_node3=>false)
        end
    end
    return isirregular
end

isirregular = refine_grid_Triangle!(grid, refineCells)