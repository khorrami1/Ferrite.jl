using Ferrite, FerriteGmsh, SparseArrays
#grid = togrid("docs/src/literate-tutorials/l.msh");
grid = generate_grid(Quadrilateral,(2,2))
grid  = ForestBWG(grid,25)
Ferrite.refine_all!(grid,1)
Ferrite.refine_all!(grid,2)

struct Elasticity
    G::Float64
    K::Float64
end

function material_routine(material::Elasticity, ε::SymmetricTensor{2})
    (; G, K) = material
    stress(ε) = 2G * dev(ε) + K * tr(ε) * one(ε)
    ∂σ∂ε, σ = gradient(stress, ε, :all)
    return σ, ∂σ∂ε
end

E = 200e3 # Young's modulus [MPa]
ν = 0.2 # Poisson's ratio [-]
material = Elasticity(E/2(1+ν), E/3(1-2ν));

function assemble_cell!(ke, fe, cellvalues, material, ue)
    fill!(ke, 0.0)
    fill!(fe, 0.0)

    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        ## For each integration point, compute strain, stress and material stiffness
        ε = function_symmetric_gradient(cellvalues, q_point, ue)
        σ, ∂σ∂ε = material_routine(material, ε)

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
            fe[i] += σ ⊡ ∇Nᵢ * dΩ
            for j in 1:n_basefuncs
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                ke[i, j] += (∂σ∂ε ⊡ ∇ˢʸᵐNⱼ) ⊡ ∇Nᵢ * dΩ
            end
        end
    end
end

function assemble_global!(K, f, a, dh, cellvalues, material)
    ## Allocate the element stiffness matrix and element force vector
    n_basefuncs = getnbasefunctions(cellvalues)
    ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    ## Create an assembler
    assembler = start_assemble(K, f)
    ## Loop over all cells
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        @views ue = a[celldofs(cell)]
        ## Compute element contribution
        assemble_cell!(ke, fe, cellvalues, material, ue)
        ## Assemble ke and fe into K and f
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    return K, f
end

function solve(grid)
    dim = 2
    order = 1 
    ip = Lagrange{RefQuadrilateral, order}()^dim 
    qr = QuadratureRule{RefQuadrilateral}(2) 
    cellvalues = CellValues(qr, ip);

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh);

    ch = ConstraintHandler(dh)
    add!(ch, Ferrite.ConformityConstraint(:u))
    # add!(ch, Dirichlet(:u, getfaceset(grid, "top"), (x, t) -> Vec{2}((0.0,0.0)), [1,2]))
    add!(ch, Dirichlet(:u, getfaceset(grid, "left"), (x, t) -> Vec{2}((-0.1,0.0)), [1,2]))
    add!(ch, Dirichlet(:u, getfaceset(grid, "right"), (x, t) -> Vec{2}((0.1,0.1)), [1,2]))
    close!(ch);

    K = create_sparsity_pattern(dh,ch)
    f = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    assemble_global!(K, f, a, dh, cellvalues, material);
    apply!(K, f, ch)
    u = K \ f;
    apply!(u,ch)
    return u,dh,ch,cellvalues
end

function compute_fluxes(u,dh)
    ip = Lagrange{RefQuadrilateral, 1}()^2
    ip_1 = Lagrange{RefQuadrilateral, 1}()
    ## Superconvergent points
    qr_1 = QuadratureRule{RefQuadrilateral}(1)
    cellvalues_sc = CellValues(qr_1, ip);
    ## "Normal" quadrature points for the fluxes
    qr_2 = QuadratureRule{RefQuadrilateral}(2)
    cellvalues = CellValues(qr_2, ip);
    ## for gradient of stress
    # cellvalues_gs = CellValues(qr_1, ip_1)
    # cellvalues_gs_sc = CellValues(qr_2, ip_1)
    ## Buffers
    σ_gp_sc = Vector{Vector{SymmetricTensor{2,2,Float64,3}}}()
    σ_gp_sc_loc = Vector{SymmetricTensor{2,2,Float64,3}}()
    σ_gp = Vector{Vector{SymmetricTensor{2,2,Float64,3}}}()
    σ_gp_loc = Vector{SymmetricTensor{2,2,Float64,3}}()
    ε_gp_sc = Vector{Vector{SymmetricTensor{2,2,Float64,3}}}()
    ε_gp_sc_loc = Vector{SymmetricTensor{2,2,Float64,3}}()
    ε_gp = Vector{Vector{SymmetricTensor{2,2,Float64,3}}}()
    ε_gp_loc = Vector{SymmetricTensor{2,2,Float64,3}}()
    # ∇σ_gp = Vector{Vector{Vec{2, Float64}}}()
    # ∇σ_gp_loc = Vector{Vec{2, Float64}}()
    # ∇σ_gp_sc = Vector{Vector{Vec{2, Float64}}}()
    # ∇σ_gp_sc_loc = Vector{Vec{2, Float64}}()
    for (cellid,cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues, cell)
        reinit!(cellvalues_sc, cell)
        @views ue = u[celldofs(cell)]
        for q_point in 1:getnquadpoints(cellvalues)
            ε = function_symmetric_gradient(cellvalues, q_point, ue)
            σ, _ = material_routine(material, ε)
            push!(σ_gp_loc, σ)
            push!(ε_gp_loc, ε)
        end
        for q_point in 1:getnquadpoints(cellvalues_sc)
            ε = function_symmetric_gradient(cellvalues_sc, q_point, ue)
            σ, _ = material_routine(material, ε)
            push!(σ_gp_sc_loc, σ)
            push!(ε_gp_sc_loc, ε)
        end
        push!(σ_gp,copy(σ_gp_loc))
        push!(σ_gp_sc,copy(σ_gp_sc_loc))
        push!(ε_gp,copy(ε_gp_loc))
        push!(ε_gp_sc,copy(ε_gp_sc_loc))
        ## Reset buffer for local points
        empty!(σ_gp_loc)
        empty!(σ_gp_sc_loc)
        empty!(ε_gp_loc)
        empty!(ε_gp_sc_loc)
    end

    # # projector = L2Projector(Lagrange{RefQuadrilateral, 1}()^2, transfered_grid)
    # projector = L2Projector(Lagrange{RefQuadrilateral, 1}()^2, dh.grid)
    # σ_dof = project(projector, σ_gp, QuadratureRule{RefQuadrilateral}(2))

    # for (cellid, cell) in enumerate(CellIterator(dh))
    #     reinit!(cellvalues, cell)
    #     reinit!(cellvalues_sc, cell)
    #     reinit!(cellvalues_gs, cell)
    #     reinit!(cellvalues_gs_sc, cell)
    #     @views σ_dof_e = σ_dof[cell.nodes]
    #     σ_dof_e_11 = map(x->x[1,1], σ_dof_e)
    #     for q_point in 1:getnquadpoints(cellvalues_gs)
    #         ∇σ = function_gradient(cellvalues_gs, q_point, σ_dof_e_11)
    #         push!(∇σ_gp_loc, ∇σ)
    #     end
    #     for q_point in 1:getnquadpoints(cellvalues_gs_sc)
    #         ∇σ = function_gradient(cellvalues_gs_sc, q_point, σ_dof_e_11)
    #         push!(∇σ_gp_sc_loc, ∇σ)
    #     end
    #     push!(∇σ_gp, copy(∇σ_gp_loc))
    #     push!(∇σ_gp_sc, copy(∇σ_gp_sc_loc))
    #     empty!(∇σ_gp_loc)
    #     empty!(∇σ_gp_sc_loc)
    # end

    return σ_gp, σ_gp_sc, ε_gp, ε_gp_sc
end

function solve_adaptive(initial_grid)
    ip = Lagrange{RefQuadrilateral, 1}()
    qr = QuadratureRule{RefQuadrilateral}(1)
    cellvalues_tensorial = CellValues(qr, ip);
    finished = false
    i = 1
    grid = initial_grid
    while !finished && i<=10
        transfered_grid = Ferrite.creategrid(grid)
        u,dh,ch,cv = solve(transfered_grid)
        σ_gp, σ_gp_sc, ε_gp, ε_gp_sc = compute_fluxes(u,dh)
        projector = L2Projector(Lagrange{RefQuadrilateral, 1}()^2, transfered_grid)
        σ_dof = project(projector, σ_gp, QuadratureRule{RefQuadrilateral}(2))
        ε_dof = project(projector, ε_gp, QuadratureRule{RefQuadrilateral}(2))
        # ∇σ_dof = project(projector, ∇σ_gp, QuadratureRule{RefQuadrilateral}(1))
        cells_to_refine = Int[]
        error_arr = Float64[]
        for (cellid,cell) in enumerate(CellIterator(projector.dh))
            reinit!(cellvalues_tensorial, cell)
            # @views σe = σ_dof[celldofs(cell)]
            @views εe = ε_dof[celldofs(cell)]
            error = 0.0
            for q_point in 1:getnquadpoints(cellvalues_tensorial)
                # σ_dof_at_sc = function_value(cellvalues_tensorial, q_point, σe)
                ε_dof_at_sc = function_value(cellvalues_tensorial, q_point, εe)
                # error += norm((σ_gp_sc[cellid][1] - σ_dof_at_sc )) * getdetJdV(cellvalues_tensorial,q_point)
                error += norm((ε_gp_sc[cellid][1] - ε_dof_at_sc )) * getdetJdV(cellvalues_tensorial,q_point)
                # error += norm((∇σ_gp_sc[cellid][1] - ∇σ_gp_loc )) * getdetJdV(cellvalues_tensorial,q_point)
            end
            if error > 1.0e-4
                push!(cells_to_refine,cellid)
            end
            push!(error_arr,error)
        end
        vtk_grid("Results/linear_elasticity-$i", dh) do vtk
            vtk_point_data(vtk, dh, u)
            vtk_point_data(vtk, projector, σ_dof, "stress")
            vtk_cell_data(vtk, getindex.(collect(Iterators.flatten(σ_gp_sc)),1), "stress sc")
            vtk_cell_data(vtk, getindex.(collect(Iterators.flatten(ε_gp_sc)),1), "strain sc")
            vtk_cell_data(vtk, error_arr, "error")
        end

        Ferrite.refine!(grid,cells_to_refine)
        vtk_grid("unbalanced.vtu", dh) do vtk
        end

        Ferrite.balanceforest!(grid)
        vtk_grid("balanced.vtu", dh) do vtk
        end

        i += 1
        if isempty(cells_to_refine)
            finished = true
        end
    end
end

solve_adaptive(grid)
