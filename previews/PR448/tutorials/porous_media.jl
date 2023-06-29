using Ferrite, FerriteMeshParser, Tensors

struct Elastic{T}
    C::SymmetricTensor{4,2,T,9}
end
function Elastic(;E=20.e3, ν=0.3)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    I2 = one(SymmetricTensor{2,2})
    I4vol = I2⊗I2
    I4dev = minorsymmetric(otimesu(I2,I2)) - I4vol / 3
    return Elastic(2G*I4dev + K*I4vol)
end

function element_routine!(Ke, fext, material::Elastic, cv, cell, args...)
    reinit!(cv, cell)
    n_basefuncs = getnbasefunctions(cv)
    dσdϵ = material.C
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            δ∇N = shape_symmetric_gradient(cv, q_point, i)
            for j in 1:n_basefuncs
                ∇N = shape_symmetric_gradient(cv, q_point, j)
                Ke[i, j] += δ∇N ⊡ dσdϵ ⊡ ∇N * dΩ
            end
        end
    end
end

struct PoroElastic{T}
    elastic::Elastic{T} ## Skeleton stiffness
    k::T            ## Permeability of liquid   [mm^4/(Ns)]
    ϕ::T            ## Porosity                 [-]
    K_liquid::T     ## Liquid bulk modulus      [MPa]
end
PoroElastic(;elastic, k, ϕ, K_liquid) = PoroElastic(elastic, k, ϕ, K_liquid) # Keyword constructor

function element_routine!(Ke, fext, m::PoroElastic, cvs::Tuple, cell, a_old, Δt, sdh)
    # Setup cellvalues and give easier names
    reinit!.(cvs, (cell,))
    cv_u, cv_p = cvs

    # Check that cellvalues are compatible with each other (should have same quadrature rule)
    @assert getnquadpoints(cv_u) == getnquadpoints(cv_p)

    C_el = m.elastic.C ## Elastic stiffness

    # Assemble stiffness and force vectors
    for q_point in 1:getnquadpoints(cv_u)
        dΩ = getdetJdV(cv_u, q_point)
        # Variation of u_i
        for (iᵤ, Iᵤ) in pairs(dof_range(sdh, :u))
            ∇δNu = shape_symmetric_gradient(cv_u, q_point, iᵤ)
            div_δNu = shape_divergence(cv_u, q_point, iᵤ)
            for (jᵤ, Jᵤ) in pairs(dof_range(sdh, :u))
                ∇Nu = shape_symmetric_gradient(cv_u, q_point, jᵤ)
                Ke[Iᵤ, Jᵤ] -= ∇δNu ⊡ C_el ⊡ ∇Nu * dΩ
            end
            for (jₚ, Jₚ) in pairs(dof_range(sdh, :p))
                Np = shape_value(cv_p, q_point, jₚ)
                Ke[Iᵤ, Jₚ] += div_δNu * Np
            end
        end
        # Variation of p_i
        for (iₚ, Iₚ) in pairs(dof_range(sdh, :p))
            δNp = shape_value(cv_p, q_point, iₚ)
            ∇δNp = shape_gradient(cv_p, q_point, iₚ)
            for (jᵤ, Jᵤ) in pairs(dof_range(sdh, :u))
                div_Nu = shape_divergence(cv_u, q_point, jᵤ)
                Lpu_ij = δNp*div_Nu*dΩ
                Ke[Iₚ,Jᵤ] += Lpu_ij
                fext[Iₚ] += Lpu_ij*a_old[Jᵤ]
            end
            for (jₚ, Jₚ) in pairs(dof_range(sdh, :p))
                ∇Np = shape_gradient(cv_p, q_point, jₚ)
                Np = shape_value(cv_p, q_point, jₚ)
                Kpp_ij = (m.k/m.ϕ) * ∇δNp ⋅ ∇Np * dΩ
                Lpp_ij = δNp*Np/m.K_liquid
                Ke[Iₚ,Jₚ] += Δt*Kpp_ij + Lpp_ij
                fext[Iₚ] += Lpp_ij*a_old[Jₚ]
            end
        end
    end
end

struct FEDomain{M,CV,SDH<:SubDofHandler}
    material::M
    cellvalues::CV
    sdh::SDH
end

function doassemble!(K, f, domains::Vector{<:FEDomain}, a_old, Δt)
    assembler = start_assemble(K, f)
    for domain in domains
        doassemble!(assembler, domain, a_old, Δt)
    end
end

function doassemble!(assembler, domain::FEDomain, a_old, Δt)
    material = domain.material
    cv = domain.cellvalues
    sdh = domain.sdh
    n = ndofs_per_cell(sdh)
    Ke = zeros(n,n)
    fe = zeros(n)
    ae_old = zeros(n)

    for cell in CellIterator(sdh)
        map!(i->a_old[i], ae_old, celldofs(cell)) # copy values from a_old to ae_old
        fill!(Ke, 0)
        fill!(fe, 0)
        element_routine!(Ke, fe, material, cv, cell, ae_old, Δt, sdh)
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
end

function get_grid()
    # Import grid from abaqus mesh
    grid = get_ferrite_grid(joinpath(@__DIR__, "porous_media_0p25.inp"))

    # Create cellsets for each fieldhandler
    addcellset!(grid, "solid3", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS3")))
    addcellset!(grid, "solid4", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS4R")))
    addcellset!(grid, "porous3", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS3")))
    addcellset!(grid, "porous4", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS4R")))
    return grid
end

function setup_problem(;t_rise=0.1, p_max=100.0)

    grid = get_grid()

    # Setup the materials
    m_solid = Elastic(;E=20.e3, ν=0.3)
    m_porous = PoroElastic(;elastic=Elastic(;E=10e3, ν=0.3), K_liquid=15e3, k=5.0e-2, ϕ=0.8)

    # Setup the interpolation and integration rules
    dim=Ferrite.getdim(grid)
    ipu_quad = Lagrange{RefQuadrilateral,2}()^2
    ipu_tri  = Lagrange{RefTriangle,2}()^2
    ipp_quad = Lagrange{RefQuadrilateral,1}()
    ipp_tri  = Lagrange{RefTriangle,1}()

    qr_quad = QuadratureRule{RefQuadrilateral}(2)
    qr_tri  = QuadratureRule{RefTriangle}(2)

    # CellValues
    cvu_quad = CellValues(qr_quad, ipu_quad)
    cvu_tri = CellValues(qr_tri, ipu_tri)
    cvp_quad = CellValues(qr_quad, ipp_quad)
    cvp_tri = CellValues(qr_tri, ipp_tri)

    # Setup the DofHandler
    dh = DofHandler(grid)
    # Solid quads
    sdh_solid_quad = SubDofHandler(dh, getcellset(grid,"solid4"))
    add!(sdh_solid_quad, :u, ipu_quad)
    # Solid triangles
    sdh_solid_tri = SubDofHandler(dh, getcellset(grid,"solid3"))
    add!(sdh_solid_tri, :u, ipu_tri)
    # Porous quads
    sdh_porous_quad = SubDofHandler(dh, getcellset(grid, "porous4"))
    add!(sdh_porous_quad, :u, ipu_quad)
    add!(sdh_porous_quad, :p, ipp_quad)
    # Porous triangles
    sdh_porous_tri = SubDofHandler(dh, getcellset(grid, "porous3"))
    add!(sdh_porous_tri, :u, ipu_tri)
    add!(sdh_porous_tri, :p, ipp_tri)

    close!(dh)

    # Setup the domains
    domains = [FEDomain(m_solid, cvu_quad, sdh_solid_quad),
               FEDomain(m_solid, cvu_tri, sdh_solid_tri),
               FEDomain(m_porous, (cvu_quad, cvp_quad), sdh_porous_quad),
               FEDomain(m_porous, (cvu_tri, cvp_tri), sdh_porous_tri)
               ]

    # Add boundary conditions
    ch = ConstraintHandler(dh);
    add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), (x, t) -> zero(Vec{2}), [1,2]))
    add!(ch, Dirichlet(:p, getfaceset(grid, "bottom_p"), (x, t) -> 0.0))
    add!(ch, Dirichlet(:p, getfaceset(grid, "top_p"), (x, t) -> p_max*clamp(t/t_rise,0,1)))
    close!(ch)

    return dh, ch, domains
end

function solve(dh, ch, domains; Δt=0.025, t_total=1.0)
    K = create_sparsity_pattern(dh);
    f = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    pvd = paraview_collection("porous_media.pvd");
    for (step, t) = enumerate(0:Δt:t_total)
        if t>0
            doassemble!(K, f, domains, a, Δt)
            update!(ch, t)
            apply!(K, f, ch)
            a .= K\f
        end
        vtk_grid("porous_media-$step", dh) do vtk
            vtk_point_data(vtk, dh, a)
            vtk_save(vtk)
            pvd[step] = vtk
        end
    end
    vtk_save(pvd);
end

dh, ch, domains = setup_problem()
solve(dh, ch, domains)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

