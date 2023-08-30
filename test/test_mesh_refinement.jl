

using Ferrite, Tensors

grid = generate_grid(Triangle, (10,10))

# refineCell = Vector{Int}
refineCells = [1, 2]

# Triangle element refinement

function refine_grid_Triangle!(grid::Grid, refineCells::Vector{Int64})

    newNodes = copy(grid.nodes)
    newCells = copy(grid.cells)

    # isirregular = Pair{Int64, Bool}[]
    isirregular = Bool[]
    addIdNodes = Int[]
    addIdNodesVicinity = Vector{Int}[]
    addNodes = Node[]

    for cellId in refineCells

        nodes_cell = grid.cells[cellId].nodes

        cellCoord = Ferrite.getcoordinates(grid, cellId)

        p1 = (cellCoord[1]+cellCoord[2])/2
        p2 = (cellCoord[1]+cellCoord[3])/2
        p3 = (cellCoord[2]+cellCoord[3])/2

        node1 = Node((p1[1], p1[2]))
        node2 = Node((p2[1], p2[2]))
        node3 = Node((p3[1], p3[2]))

        nnodes = length(newNodes) #number of total nodes including the new nodes

        node1_temp = findall(x->x== node1, addNodes)

        if isempty(node1_temp)
            id_node1 = nnodes + 1
            push!(addIdNodes, id_node1)
            push!(addIdNodesVicinity, [nodes_cell[1], nodes_cell[2]])
            push!(newNodes, node1)
            push!(isirregular, true)
            push!(addNodes, node1)
        else
            nnodes -= 1
            idNode_temp = node1_temp[1]
            isirregular[idNode_temp] = false
            id_node1 = addIdNodes[idNode_temp]
        end

        node2_temp = findall(x->x== node2, addNodes)

        if isempty(node2_temp)
            id_node2 = nnodes + 2
            push!(addIdNodes, id_node2)
            push!(addIdNodesVicinity, [nodes_cell[1], nodes_cell[3]])
            push!(newNodes, node2)
            push!(isirregular, true)
            push!(addNodes, node2)
        else
            nnodes -= 1
            idNode_temp = node2_temp[1]
            isirregular[idNode_temp] = false
            id_node2 = addIdNodes[idNode_temp]
        end

        node3_temp = findall(x->x== node3, addNodes)

        if isempty(node3_temp)
            id_node3 = nnodes + 3
            push!(addIdNodes, id_node3)
            push!(addIdNodesVicinity, [nodes_cell[2], nodes_cell[3]])
            push!(newNodes, node3)
            push!(isirregular, true)
            push!(addNodes, node3)
        else
            nnodes -= 1
            idNode_temp = node3_temp[1]
            isirregular[idNode_temp] = false
            id_node3 = addIdNodes[idNode_temp]
        end

        push!(newCells, Triangle((id_node1, id_node2, id_node3)))

    end

    return newNodes, newCells, isirregular, addIdNodes, addIdNodesVicinity
end

newNodes, newCells, isirregular, addIdNodes, addIdNodesVicinity = refine_grid_Triangle!(grid, refineCells)

newGrid = Grid(newCells, newNodes)
newDH = DofHandler(newGrid)
add!(newDH, :u, 1)
close!(newDH)

ch = Ferrite.ConstraintHandler(newDH)

for i in 1:length(addIdNodes)
    if isirregular[i]
        # Adding AffineConstraint for the corresponding dofs of the corresponding node
        idNode = addIdNodes[i]
        idNodeVicinity = addIdNodesVicinity[i]
        dof = idNode*1 # if dim=1 (add!(newDH, :u, 1))
        affineConstraint = AffineConstraint(dof, [idNodeVicinity[1]=>0.5, idNodeVicinity[2]=>0.5], 0.0)
        add!(ch, affineConstraint)
    end
end
