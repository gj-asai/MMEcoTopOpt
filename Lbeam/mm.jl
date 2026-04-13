"""
L beam optimization - multimaterial
CO2 minimization with compliance and total volume fraction constraints
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

function mm(comp_max, volfrac, rρ, rθ, angle=0; filename)
    reset_timer!()
    pvd = paraview_collection(filename)

    # read mesh
    mesh_file = "models/Lbeam.msh"
    grid = redirect_stdout(() -> togrid(mesh_file), devnull) # suppress print from the Gmsh.jl call inside togrid
    @info "Done reading $(mesh_file): $(getnnodes(grid)) nodes, $(getncells(grid)) elements"

    # loads and BCs
    Lsize = 50
    addfacetset!(grid, "fixed", x -> x[2] ≈ 2*Lsize) # top edge
    addnodeset!(grid, "force", x -> x[1] ≈ 2*Lsize && x[2] ≈ Lsize) # top right corner

    constraints = [
        Dirichlet(:u, getfacetset(grid, "fixed"), (x, t) -> [0.0, 0.0]), # block x and y
    ]
    loads = [
        NodalLoad("force", (0.0, -100.0)),
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
    model = FEModel(; grid, ip, qr, mat_interp, constraints)
    f = compute_rhs(loads, model)

    # initialize design variables
    x = collect(Iterators.flatten(zip(
        fill(volfrac, getncells(grid)), # ρ1
        fill(volfrac, getncells(grid)), # ρ2
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

    dCO2dx, dcdx, dvdx = similar(x), similar(x), similar(x)
    dgdx = zeros(2, 3 * getncells(grid))

    # initialize linear solver
    fesolver = FESolver(model)

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
        ρi = @view x[i:nmaterials+1:end]
        CO2_ini += mati.CO2 * mati.ρ * ρi ⋅ model.elemvol
    end

    history = Dict(
        :beta => Float64[],
        :penal => Float64[],
        :objective => Float64[],
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
                xTilde .= x
                for i in 1:nmaterials
                    @views TopOptTools.filter!(xTilde[i:nmaterials+1:end], density_filter)
                end
                @views TopOptTools.filter!(xTilde[nmaterials+1:nmaterials+1:end], orientation_filter)

                xPhys .= xTilde
                for i in 1:nmaterials
                    @views project_heaviside!(xPhys[i:nmaterials+1:end], beta, eta)
                end
            end

            @timeit "assemble stiffness" update_stiffness!(fesolver, xPhys)
            @timeit "linear solve" begin
                fea!(fesolver, f)
                u = fesolver.solution
            end

            @timeit "evaluate functions" begin
                # Objective: CO2 impact
                CO2 = 0
                for i in 1:nmaterials
                    mati = mat_interp.mat[i]
                    @views ρi = xPhys[i:nmaterials+1:end]
                    CO2 += mati.CO2 * mati.ρ * ρi ⋅ model.elemvol
                end

                # Constraint: max compliance
                c = dot(u, fesolver.K, u)

                # Constraint: max volume fraction
                v = 0
                for i in 1:nmaterials
                    @views v += xPhys[i:nmaterials+1:end] ⋅ model.elemvol
                end

                # all constraints
                g = [c / comp_max - 1, v / sum(model.elemvol) / volfrac - 1]
            end

            @timeit "sensitivity analysis" begin
                # Objective: CO2 impact
                fill!(dCO2dx, 0.0)
                for i in 1:nmaterials
                    mati = mat_interp.mat[i]
                    dCO2dx[i:nmaterials+1:end] .= mati.CO2 * mati.ρ * model.elemvol
                end

                # Constraint: max compliance
                adjoint_sensitivities!(dcdx, u, u, fesolver)

                # Constraint: max volume fraction
                fill!(dvdx, 0.0)
                for i in 1:nmaterials
                    dvdx[i:nmaterials+1:end] .= model.elemvol
                end

                # Chain rule
                dPhysdTilde .= xTilde
                @views project_heaviside_derivative!(dPhysdTilde, beta, eta)
                # densities
                for i in 1:nmaterials
                    dCO2dx[i:nmaterials+1:end] .*= dPhysdTilde[i:nmaterials+1:end] # xPhys -> xTilde
                    @views TopOptTools.filter!(dCO2dx[i:nmaterials+1:end], density_filter) # xTilde -> x

                    dcdx[i:nmaterials+1:end] .*= dPhysdTilde[i:nmaterials+1:end] # xPhys -> xTilde
                    @views TopOptTools.filter!(dcdx[i:nmaterials+1:end], density_filter) # xTilde -> x

                    dvdx[i:nmaterials+1:end] .*= dPhysdTilde[i:nmaterials+1:end] # xPhys -> xTilde
                    @views TopOptTools.filter!(dvdx[i:nmaterials+1:end], density_filter) # xTilde -> x
                end
                # angles, no need to filter dCO2dθ and dvdθ, they are all zero
                @views TopOptTools.filter!(dcdx[nmaterials+1:nmaterials+1:end], orientation_filter)

                dgdx[1, :] .= dcdx / comp_max
                dgdx[2, :] .= dvdx / sum(model.elemvol) / volfrac
            end

            # volume fractions of each material
            volfracs = Float64[]
            for i = 1:nmaterials
                push!(volfracs, xPhys[i:nmaterials+1:end] ⋅ model.elemvol / sum(model.elemvol))
            end

            # Stopping criterion: relative change in the objective function
            change = if loop > 1
                abs(CO2 - history[:objective][end]) / history[:objective][end]
            else
                Inf
            end

            formatted_volfracs = join([@sprintf("%5.2f", 100 * f) for f in volfracs], ", ")
            @info @sprintf "It = %4d | CO2 = %8.4f | change = %8.2e | c = %9.4f | volfracs = [%s] %%" loop CO2 change c formatted_volfracs

            # update history
            push!(history[:beta], beta)
            push!(history[:penal], mat_interp.penal)
            push!(history[:objective], CO2)
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

                    write_solution(vtk, model.dh, fesolver.solution)
                    pvd[loop] = vtk
                end
            end

            # apply continuation
            if (mma.iter >= 50 && change < 1e-4) || mma.iter >= 150
                # if at max beta, stop the optimization
                beta == 64 && break
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
        print_timer(title="c=$(comp_max), v=$(volfrac), angle=$(angle)")
    end
end
