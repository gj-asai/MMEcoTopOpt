"""
MBB optimization - multimaterial
CO2 minimization with displacement and mass constraints
"""

using Ferrite, FerriteGmsh, UnPack
using TimerOutputs, Printf, WriteVTK, JLD2, HDF5
using TopOptTools

# define interpolation scheme
# struct needs to be mutable because we want to update penal during the optimization
mutable struct MMSOMP{dim,M,N,T<:Real,CT} <: MaterialInterpolation{N,T}
    mat::Vector{Material{dim,T,CT}}
    penal::T
end
MMSOMP(mat::Vector{Material{dim,T,CT}}, penal::T) where {dim,T,CT} = MMSOMP{dim,length(mat),length(mat) + 1,T,CT}(mat, penal)

# UMMI-2, for an arbitrary number of materials in 2D and in 3D
function TopOptTools.interpolate(xe::AbstractVector{T}, interp::MMSOMP{dim}) where {T<:Real,dim}
    ρ, θ = xe[1:end-1], xe[end]
    result = zero(SymmetricTensor{4,2,T})
    for i in eachindex(interp.mat)
        weight = one(T)
        for j in eachindex(interp.mat)
            if i == j
                weight *= ρ[j]^interp.penal
            else
                weight *= one(T) - ρ[j]^interp.penal
            end
        end
        result += weight * interp.mat[i].C
    end
    result = TopOptTools.void(Val(dim)).C + rotate(result, θ)
    return result
end
rotate_mmsomp(mat::Material{2}, θ) = rotate(mat.C, θ)
rotate_mmsomp(mat::Material{3}, θ) = rotate(mat.C, Vec{3}((0.0, 0.0, 1.0)), θ)

function mm_alternative(disp_max, mass_max, rρ, rθ, angle=0; filename)
    # constraints are given for the full beam with 8 mm thickness
    # converting them for a half beam with unit thickness
    mass_max /= 16.0
    disp_max *= 8.0

    reset_timer!()
    pvd = paraview_collection(filename)

    # read mesh
    mesh_file = "models/mbb.msh" # this example contains a rectangular mesh
    grid = redirect_stdout(() -> togrid(mesh_file), devnull) # suppress print from the Gmsh.jl call inside togrid
    @info "Done reading $(mesh_file): $(getnnodes(grid)) nodes, $(getncells(grid)) elements"

    # loads and BCs
    xmax = maximum([n.x[1] for n in getnodes(grid)])
    ymax = maximum([n.x[2] for n in getnodes(grid)])
    addfacetset!(grid, "symmetry", x -> x[1] ≈ 0.0) # left edge
    addnodeset!(grid, "support", x -> x[1] ≈ xmax && x[2] ≈ 0.0) # bottom right corner
    addnodeset!(grid, "force", x -> x[1] ≈ 0.0 && x[2] ≈ ymax) # top left corner

    constraints = [
        Dirichlet(:u, getfacetset(grid, "symmetry"), (x, t) -> 0.0, [1]), # block x displacement
        Dirichlet(:u, getnodeset(grid, "support"), (x, t) -> 0.0, [2]), # block y displacement
    ]
    loads = [
        NodalLoad("force", (0.0, -700.0)),
    ]

    # 2 materials but the rest of the code works as-is if you add more materials here
    carbon = Orthotropic2D(El=16.7e3, Et=3.7e3, nult=0.33, Glt=1.4e3, ρ=1.27e-3, CO2=10.1)
    bamboo = Orthotropic2D(El=10.8e3, Et=4.6e3, nult=0.36, Glt=1.7e3, ρ=1.12e-3, CO2=2.9)
    mat_interp = MMSOMP([carbon, bamboo], 1.0)
    nmaterials = length(mat_interp.mat)

    # element interpolation and quadrature
    ip = Lagrange{RefQuadrilateral,1}() # linear elements
    qr = QuadratureRule{RefQuadrilateral}(2) # 2 point quadrature

    # FE model
    model = FEModel(; grid, ip, qr, constraints)
    f = compute_rhs(loads, model)

    # L ⋅ u = y displacement at node "force"
    L = zeros(ndofs(model.dh))
    disp_node = getnodeset(grid, "force")[1]
    for cell in CellIterator(model.dh), (i, node) in enumerate(getnodes(cell))
        node == disp_node || continue
        dofs = celldofs(cell)
        L[dofs[2i]] = -1 # y displacement, flipping the sign
        break
    end

    # initialize design variables
    x = collect(Iterators.flatten(zip(
        fill(0.5, getncells(grid)), # ρ1
        fill(0.5, getncells(grid)), # ρ2
        fill(deg2rad(angle), getncells(grid)) # θ
    )))
    xmin = collect(Iterators.flatten(zip(
        fill(0, getncells(grid)), # ρ1
        fill(0, getncells(grid)), # ρ2
        fill(-π, getncells(grid)) # θ
    )))
    xmax = collect(Iterators.flatten(zip(
        fill(1, getncells(grid)), # ρ1
        fill(1, getncells(grid)), # ρ2
        fill(π, getncells(grid)) # θ
    )))

    xTilde = copy(x)
    xPhys = copy(x)
    dPhysdTilde = similar(x)

    dCO2dx, dmdx, dudx = similar(x), similar(x), similar(x)
    dgdx = zeros(2, 3 * getncells(grid))

    # initialize linear solver
    fesolver = FEA(model, LinearElasticity(mat_interp))

    @timeit "build filters" begin
        @info "Building density filter with radius $rρ"
        density_filter = ConvolutionFilter(rρ, model)
        @info "Building orientation filter with radius $rθ"
        orientation_filter = ConvolutionFilter(rθ, model)
    end

    # initial impact, will be the normalization factor of the objective function
    CO2_ini = 0
    for i in 1:nmaterials
        mati = mat_interp.mat[i]
        ρi = x[i:nmaterials+1:end]
        CO2_ini += mati.CO2 * mati.ρ * ρi ⋅ model.elemvol
    end

    history = Dict(
        :beta => Float64[],
        :penal => Float64[],
        :objective => Float64[],
        :mass => Float64[],
        :disp => Float64[],
        :constraint => Vector{Float64}[],
        :volfracs => Vector{Float64}[],
    )

    # initialize MMA
    m, n = 2, length(x)
    a0mma, amma, cmma, dmma = 1.0, zeros(m), fill(1000.0, m), zeros(m)
    mma = MMAWorkspace(m, n, xmin, xmax, a0mma, amma, cmma, dmma, asyinit=0.1, asyincr=1.1, asydecr=0.6, move=0.5)

    beta = 1 # thresholding projection sharpness
    eta = 0.5 # threshold
    @info "Starting optimization with beta = $(beta) and p = $(mat_interp.penal)"
    try
        maxiter = 1200
        for loop in 1:maxiter+1
            @timeit "filter and project" begin
                # filter
                xTilde .= x
                for i in 1:nmaterials
                    @views TopOptTools.filter!(xTilde[i:nmaterials+1:end], density_filter)
                end
                @views TopOptTools.filter!(xTilde[nmaterials+1:nmaterials+1:end], orientation_filter)

                # project
                for i in 1:nmaterials
                    @. xPhys[i:nmaterials+1:end] = (tanh(beta * eta) + tanh(beta * (xTilde[1:3:end] - eta))) / (tanh(beta * eta) + tanh(beta * (1 - eta)))
                end
                @. xPhys[nmaterials+1:nmaterials+1:end] = xTilde[nmaterials+1:nmaterials+1:end]
            end

            # FE analysis
            @timeit "assemble stiffness matrix" update_stiffness!(fesolver, xPhys)
            @timeit "linear solve" solve!(fesolver, f)

            @timeit "evaluate functions" begin
                # CO2 impact and mass
                CO2 = 0
                mass = 0
                volfracs = Float64[]
                for i in 1:nmaterials
                    mati = mat_interp.mat[i]
                    @views ρi = xPhys[i:nmaterials+1:end]

                    CO2 += mati.CO2 * mati.ρ * ρi ⋅ model.elemvol
                    mass += mati.ρ * ρi ⋅ model.elemvol

                    push!(volfracs, ρi ⋅ model.elemvol / sum(model.elemvol))
                end

                # Constraint: displacement at node "force"
                u = copy(fesolver.solution)
                disp = L ⋅ u

                g = [disp / disp_max - 1, mass / mass_max - 1]
            end

            @timeit "sensitivity analysis" begin
                # CO2 impact and mass
                fill!(dCO2dx, 0.0)
                fill!(dmdx, 0.0)
                for i in 1:nmaterials
                    mati = mat_interp.mat[i]
                    @views ρi = xPhys[i:nmaterials+1:end]

                    dCO2dx[i:nmaterials+1:end] .= mati.CO2 * mati.ρ * model.elemvol
                    dmdx[i:nmaterials+1:end] .= model.elemvol * mat_interp.mat[i].ρ
                end

                # Constraint: displacement at node "force"
                solve!(fesolver, L)
                lambda = fesolver.solution
                adjoint_sensitivities!(dudx, lambda, u, fesolver)

                # Chain rule - xPhys -> xTilde
                for i in 1:nmaterials
                    @. dPhysdTilde[i:nmaterials+1:end] = beta * sech(beta * (xTilde[i:nmaterials+1:end] - eta))^2 / (tanh(beta * eta) + tanh(beta * (1 - eta)))
                end
                @. dPhysdTilde[nmaterials+1:nmaterials+1:end] = 1.0
                dCO2dx .*= dPhysdTilde
                dudx .*= dPhysdTilde
                dmdx .*= dPhysdTilde

                # Chain rule - xTilde -> x
                for i in 1:nmaterials
                    @views TopOptTools.filter!(dCO2dx[i:nmaterials+1:end], density_filter)
                    @views TopOptTools.filter!(dudx[i:nmaterials+1:end], density_filter)
                    @views TopOptTools.filter!(dmdx[i:nmaterials+1:end], density_filter)
                end
                # angles, no need to filter dCO2dθ and dmdθ, they are all zero
                @views TopOptTools.filter!(dudx[nmaterials+1:nmaterials+1:end], orientation_filter)

                dgdx[1, :] .= dudx / disp_max
                dgdx[2, :] .= dmdx / mass_max
            end

            # Stopping criterion: relative change in the objective function
            change = if loop > 1
                abs(CO2 - history[:objective][end]) / history[:objective][end]
            else
                Inf
            end

            formatted_volfracs = join([@sprintf("%5.2f", 100 * f) for f in volfracs], ", ")
            @info @sprintf "It = %4d | CO2 = %8.4f | disp = %8.4f | m = %8.4f | volfracs = [%s] %% | change = %8.2e" loop CO2 disp / 8 16 * mass formatted_volfracs change

            # update history
            push!(history[:beta], beta)
            push!(history[:penal], mat_interp.penal)
            push!(history[:objective], CO2)
            push!(history[:mass], mass)
            push!(history[:disp], disp)
            push!(history[:constraint], g)
            push!(history[:volfracs], volfracs)

            # export iteration file
            @timeit "export" begin
                filename_i = @sprintf "%s.%4.4d.vtu" filename (loop - 1)
                VTKGridFile(filename_i, grid) do vtk
                    for i in 1:nmaterials
                        write_cell_data(vtk, x[i:nmaterials+1:end], "x$(i)")
                    end
                    write_cell_data(vtk, x[nmaterials+1:nmaterials+1:end], "theta")

                    for i in 1:nmaterials
                        write_cell_data(vtk, xTilde[i:nmaterials+1:end], "xTilde$(i)")
                    end
                    write_cell_data(vtk, xTilde[nmaterials+1:nmaterials+1:end], "thetaTilde")

                    for i in 1:nmaterials
                        write_cell_data(vtk, xPhys[i:nmaterials+1:end], "xPhys$(i)")
                    end
                    write_cell_data(vtk, xPhys[nmaterials+1:nmaterials+1:end], "thetaPhys")

                    write_solution(vtk, model.dh, u)
                    pvd[loop] = vtk
                end
            end

            # apply continuation
            if (mma.iter >= 50 && change < 1e-4) || mma.iter >= 150
                # if at max beta, stop the optimization
                beta == 32 && break
                # else, increase beta and p
                beta *= 2
                mat_interp.penal = min(3.0, mat_interp.penal + 1.0)
                @info "Updated beta to $(beta) and p to $(mat_interp.penal)"
                restart!(mma)
            end

            # MMA update
            @timeit "mma update" x .= mma_update!(mma, x, dCO2dx / CO2_ini, g, dgdx)
        end
    catch e
        @warn "Computation interrupted - $(typeof(e))"
        # rethrow() # uncomment this to get the stack trace for debugging
    finally
        @timeit "export" begin
            # history
            h5open(filename * ".h5", "w") do file
                write(file, "beta", history[:beta])
                write(file, "penal", history[:penal])
                write(file, "objective", history[:objective])
                write(file, "disp", history[:disp])
                write(file, "mass", history[:mass])
                write(file, "constraint", reduce(hcat, history[:constraint]))
                write(file, "volfracs", reduce(hcat, history[:volfracs]))
            end
            @info "Saved file $(filename).h5"

            # x
            save("$(filename).x.jld2", "x", x)
            @info "Saved file $(filename).x.jld2"

            # xTilde
            save("$(filename).xTilde.jld2", "xTilde", xTilde)
            @info "Saved file $(filename).xTilde.jld2"

            # xPhys
            save("$(filename).xPhys.jld2", "xPhys", xPhys)
            @info "Saved file $(filename).xPhys.jld2"

            # paraview .pvd
            vtk_save(pvd)
            @info "Saved file $(filename).pvd"
        end

        # timing
        print_timer(title="disp=$(disp_max), m=$(mass_max), angle=$(angle)")
    end
end
