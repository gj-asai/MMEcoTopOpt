"""
Bridge optimization - single material
Compliance minimization with volume fraction constraint
"""

using Ferrite, FerriteGmsh, UnPack
using TimerOutputs, Printf, WriteVTK, JLD2, HDF5
using TopOptTools

# define interpolation scheme
# struct needs to be mutable because we want to update penal during the optimization
mutable struct SOMP{dim,T<:Real,CT} <: MaterialInterpolation{2,T}
    mat::Material{dim,T,CT}
    penal::T
end
TopOptTools.interpolate(xe::AbstractVector, somp::SOMP) = TopOptTools.void(Val(2)).C + xe[1]^somp.penal * rotate(somp.mat.C, xe[2])

# main function
function somp(volfrac, rρ, rθ, angle=0; filename)
    reset_timer!()
    pvd = paraview_collection(filename)

    # read mesh
    mesh_file = "models/bridge.msh" # this example contains a rectangular mesh
    grid = redirect_stdout(() -> togrid(mesh_file), devnull) # suppress print from the Gmsh.jl call inside togrid
    @info "Done reading $(mesh_file): $(getnnodes(grid)) nodes, $(getncells(grid)) elements"

    # loads and BCs
    width = maximum([n.x[1] for n in getnodes(grid)])
    height = maximum([n.x[2] for n in getnodes(grid)])
    addfacetset!(grid, "symmetry", x -> x[1] ≈ 0.0) # left edge
    addfacetset!(grid, "force", x -> x[2] ≈ height) # top edge
    addnodeset!(grid, "fixed", x -> x[1] ≈ width && x[2] ≈ 0.0) # bottom right corner

    constraints = [
        Dirichlet(:u, getfacetset(grid, "symmetry"), (x, t) -> 0.0, [1]), # block x displacement
        Dirichlet(:u, getnodeset(grid, "fixed"), (x, t) -> [0.0, 0.0]), # block x and y
    ]
    loads = [
        LinearLoad("force", (0.0, -10.0)),
    ]

    # materials
    carbon = Orthotropic2D(El=16.7e3, Et=3.7e3, nult=0.33, Glt=1.4e3, ρ=1.27e-3, CO2=10.1)
    bamboo = Orthotropic2D(El=10.8e3, Et=4.6e3, nult=0.36, Glt=1.7e3, ρ=1.12e-3, CO2=2.9)
    mat_interp = SOMP(carbon, 1.0)

    # element interpolation and quadrature
    ip = Lagrange{RefQuadrilateral,1}() # linear elements
    qr = QuadratureRule{RefQuadrilateral}(2) # 2 point quadrature

    # FE model
    model = FEModel(; grid, ip, qr, mat_interp, constraints)
    f = compute_rhs(loads, model)

    # initialize design variables
    x = collect(Iterators.flatten(zip(
        fill(volfrac, getncells(grid)), # ρ
        fill(deg2rad(angle), getncells(grid)) # θ
    )))
    xmin = collect(Iterators.flatten(zip(
        fill(0, getncells(grid)), # ρ
        fill(-π, getncells(grid)) # θ
    )))
    xmax = collect(Iterators.flatten(zip(
        fill(1, getncells(grid)), # ρ
        fill(π, getncells(grid)) # θ
    )))

    xTilde = copy(x)
    xPhys = copy(x)
    dPhysdTilde = similar(x)

    dcdx, dvdx = similar(x), similar(x)
    dgdx = zeros(1, 2 * getncells(grid))

    # passive region
    for i in 1:getncells(grid)
        if model.centers[2, i] >= height - 3.0
            x[2i-1] = 1
            xmin[2i-1] = 0.99
        end
    end

    # initialize linear solver
    fesolver = FESolver(model)

    @timeit "build filters" begin
        @info "Building density filter"
        density_filter = ConvolutionFilter(rρ, model)
        @info "Building orientation filter"
        orientation_filter = ConvolutionFilter(rθ, model)
    end

    history = Dict(
        :beta => Float64[],
        :penal => Float64[],
        :objective => Float64[],
        :constraint => Vector{Float64}[],
        :CO2 => Float64[],
    )

    # initialize MMA
    m, n = 1, length(x)
    a0mma, amma, cmma, dmma = 1.0, zeros(m), fill(1000.0, m), zeros(m)
    mma = MMAWorkspace(m, n, xmin, xmax, a0mma, amma, cmma, dmma, asyinit=0.1, asyincr=1.1, asydecr=0.6, move=0.5)

    beta = 1
    eta = 0.5
    @info "Starting optimization with beta = $(beta) and p = $(mat_interp.penal)"
    try
        maxiter = 1000
        for loop in 1:maxiter+1
            @timeit "filter and project" begin
                xTilde .= x
                @views TopOptTools.filter!(xTilde[1:2:end], density_filter)
                @views TopOptTools.filter!(xTilde[2:2:end], orientation_filter)

                xPhys .= xTilde
                @views project_heaviside!(xPhys[1:2:end], beta, eta)
            end

            @timeit "assemble stiffness" update_stiffness!(fesolver, xPhys)
            @timeit "linear solve" begin
                fea!(fesolver, f)
                u = fesolver.solution
            end

            @timeit "evaluate functions" begin
                CO2 = mat_interp.mat.CO2 * mat_interp.mat.ρ * xPhys[1:2:end] ⋅ model.elemvol

                # Objective: max compliance
                c = dot(u, fesolver.K, u)

                # Constraint: max volume fraction
                @views v = xPhys[1:2:end] ⋅ model.elemvol

                # all constraints
                g = [v / sum(model.elemvol) / volfrac - 1,]
            end

            @timeit "sensitivity analysis" begin
                # Objective: max compliance
                adjoint_sensitivities!(dcdx, u, u, fesolver)

                # Constraint: max volume fraction
                fill!(dvdx, 0.0)
                dvdx[1:2:end] .= model.elemvol

                # Chain rule
                dPhysdTilde .= xTilde
                @views project_heaviside_derivative!(dPhysdTilde, beta, eta)

                # densities
                dcdx[1:2:end] .*= dPhysdTilde[1:2:end] # xPhys -> xTilde
                @views TopOptTools.filter!(dcdx[1:2:end], density_filter) # xTilde -> x

                dvdx[1:2:end] .*= dPhysdTilde[1:2:end] # xPhys -> xTilde
                @views TopOptTools.filter!(dvdx[1:2:end], density_filter) # xTilde -> x

                # angles, no need to filter dvdtheta, it is all zero
                @views TopOptTools.filter!(dcdx[2:2:end], orientation_filter)

                dgdx[1, :] .= dvdx / sum(model.elemvol) / volfrac
            end

            # Stopping criterion: relative change in the objective function
            change = if loop > 1
                abs(c - history[:objective][end]) / history[:objective][end]
            else
                Inf
            end

            @info @sprintf "It = %4d | CO2 = %8.4f | change = %8.2e | c = %9.4f" loop CO2 change c

            # update history
            push!(history[:beta], beta)
            push!(history[:penal], mat_interp.penal)
            push!(history[:objective], c)
            push!(history[:constraint], g)
            push!(history[:CO2], CO2)

            # export iteration file
            @timeit "export" begin
                filename_i = @sprintf "%s.%4.4d.vtu" filename (loop - 1)
                VTKGridFile(filename_i, grid) do vtk
                    write_cell_data(vtk, xPhys[1:2:end], "density")
                    write_cell_data(vtk, xPhys[2:2:end], "theta")
                    write_solution(vtk, model.dh, u)
                    pvd[loop] = vtk
                end
            end

            # apply continuation
            if (mma.iter >= 50 && change < 1e-4) || mma.iter >= 150
                # if at max beta, stop the optimization
                beta == 32 && break
                # else, increase beta and p, decrease gamma
                beta *= 2
                mat_interp.penal = min(3.0, mat_interp.penal + 1.0)
                @info "Updated beta to $(beta) and p to $(mat_interp.penal)"
                restart!(mma)
            end

            # MMA update
            @timeit "mma update" x .= mma_update!(mma, x, dcdx, g, dgdx)
        end
    catch e
        @warn "Computation interrupted - $(typeof(e))"
        # rethrow()
    finally
        @timeit "export" begin
            # history
            h5open(filename * ".h5", "w") do file
                write(file, "beta", history[:beta])
                write(file, "penal", history[:penal])
                write(file, "objective", history[:objective])
                write(file, "constraint", reduce(hcat, history[:constraint]))
                write(file, "CO2", reduce(hcat, history[:CO2]))
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
        print_timer(title="v=$(volfrac), angle=$(angle)")
    end
end
