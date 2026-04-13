"""
reads the final state from the results and re-runs the FEA to save more detailed data
called as post_fea_mm("result.xPhys.jld2")
"""

using JLD2
using Tensors, Ferrite, FerriteGmsh, UnPack
using TopOptTools

function post_fea_mm(filename)
    xPhys = load(filename)["xPhys"]

    # read mesh
    mesh_file = "models/mbb.msh" # this example contains a rectangular mesh
    grid = redirect_stdout(() -> togrid(mesh_file), devnull) # suppress print from the Gmsh.jl call inside togrid

    # loads and BCs
    xmax = maximum([n.x[1] for n in getnodes(grid)])
    ymax = maximum([n.x[2] for n in getnodes(grid)])
    addfacetset!(grid, "symmetry", x -> x[1] ≈ 0.0) # left edge
    addnodeset!(grid, "support", x -> x[1] ≈ xmax - 5.0 && x[2] ≈ 0.0) # bottom right corner
    addnodeset!(grid, "force", x -> x[1] ≈ 0.0 && x[2] ≈ ymax) # top left corner

    constraints = [
        Dirichlet(:u, getfacetset(grid, "symmetry"), (x, t) -> 0.0, [1]), # block x displacement
        Dirichlet(:u, getnodeset(grid, "support"), (x, t) -> 0.0, [2]), # block y displacement
    ]
    loads = [
        NodalLoad("force", (0.0, -700.0 / 8.0)),
    ]

    # materials
    carbon = Orthotropic2D(El=16.7e3, Et=3.7e3, nult=0.33, Glt=1.4e3, ρ=1.27e-3, CO2=10.1)
    bamboo = Orthotropic2D(El=10.8e3, Et=4.6e3, nult=0.36, Glt=1.7e3, ρ=1.12e-3, CO2=2.9)
    mat_interp = MMSOMP([carbon, bamboo], 3.0)

    # element interpolation and quadrature
    ip = Lagrange{RefQuadrilateral,1}() # linear elements
    qr = QuadratureRule{RefQuadrilateral}(2) # 2 point quadrature

    # FE model
    model = FEModel(; grid, ip, qr, mat_interp, constraints)
    f = compute_rhs(loads, model)

    # initialize linear solver
    fesolver = FESolver(model)

    # FE analysis
    update_stiffness!(fesolver, xPhys)
    fea!(fesolver, f)

    qp_strain, qp_principal, qp_principaldir = get_strains(fesolver)
    projector = L2Projector(ip^2, grid)

    new_filename = join(split(filename, '/')[1:end-1], '/') * "/fea.vtu"
    VTKGridFile(new_filename, grid) do vtk
        write_cell_data(vtk, xPhys[1:3:end], "density1")
        write_cell_data(vtk, xPhys[2:3:end], "density2")
        write_cell_data(vtk, xPhys[3:3:end], "theta")
        write_solution(vtk, model.dh, fesolver.u)

        write_projection(vtk, projector, project(projector, qp_strain, qr), "strain")
        write_projection(vtk, projector, project(projector, qp_principal, qr), "principal strain")

        for i = 1:2
            principaldir_i = [[Tensor{1,2}(vectors[:, i]) for vectors in el] for el in qp_principaldir]
            write_projection(vtk, projector, project(projector, principaldir_i, qr), "principal strain direction $i")
        end
    end
    @info "Saved $(new_filename)"
end

# calculate strains at all quadrature points
function get_strains(solver::FESolver)
    @unpack cellvalues, grid, dh = solver.model
    dim = TopOptTools.get_dim(solver.model)

    qp_global_strain = [
        [zero(SymmetricTensor{2,dim}) for _ in 1:getnquadpoints(cellvalues)]
        for _ in 1:getncells(grid)
    ]
    qp_principal = [
        [zero(SymmetricTensor{2,dim}) for _ in 1:getnquadpoints(cellvalues)]
        for _ in 1:getncells(grid)
    ]
    qp_directions = [
        [zeros(dim, dim) for _ in 1:getnquadpoints(cellvalues)]
        for _ in 1:getncells(grid)
    ]

    for cell in CellIterator(dh)
        Ferrite.reinit!(cellvalues, cell)
        e = cellid(cell)
        cell_global_strain = qp_global_strain[e]
        cell_principal = qp_principal[e]
        cell_directions = qp_directions[e]
        for q_point in 1:getnquadpoints(cellvalues)
            ε = function_symmetric_gradient(cellvalues, q_point, solver.u, celldofs(cell))
            cell_global_strain[q_point] = ε

            # sort eigenvalues by absolute value, largest first
            principal = eigen(ε)
            principal_order = sortperm(abs.(principal.values), rev=true)
            cell_principal[q_point] = diagm(SymmetricTensor{2,dim}, principal.values[principal_order])
            cell_directions[q_point] = principal.vectors[:, principal_order]
        end
    end

    return qp_global_strain, qp_principal, qp_directions
end

# this is the same as in mm.jl
mutable struct MMSOMP{dim,M,N,T<:Real,CT} <: MaterialInterpolation{N,T}
    mat::Vector{Material{dim,T,CT}}
    penal::T
end
MMSOMP(mat::Vector{Material{dim,T,CT}}, penal::T) where {dim,T,CT} = MMSOMP{dim,length(mat),length(mat) + 1,T,CT}(mat, penal)

function TopOptTools.interpolate(xe::AbstractVector{T}, interp::MMSOMP{dim}) where {T<:Real,dim}
    ρ, θ = xe[1:end-1], xe[end]
    result = zero(SymmetricTensor{4,dim,T})
    for i in eachindex(interp.mat)
        weight = one(T)
        for j in eachindex(interp.mat)
            if i == j
                weight *= ρ[j]^interp.penal
            else
                weight *= one(T) - ρ[j]^interp.penal
            end
        end
        result += weight * rotate_mmsomp(interp.mat[i], θ)
    end
    result += TopOptTools.void(Val(dim)).C
    return result
end
rotate_mmsomp(mat::Material{2}, θ) = rotate(mat.C, θ)
rotate_mmsomp(mat::Material{3}, θ) = rotate(mat.C, Vec{3}((0.0, 0.0, 1.0)), θ)
