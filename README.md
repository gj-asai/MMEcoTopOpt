# Multimaterial topology optimization of 3D printed structures

Scripts for ecoefficient optimization of multimaterial structures.
The main dependency is the package [`TopOptTools`](https://github.com/gj-asai/TopOptTools.jl.git), which is included in this repository as a Git submodule.

## Structure
The root of this repository contains the scripts for the optimization of a multimaterial MBB beam:
-  `sm.jl`: single material optimization. It is used to define the feasible limits (compliance and emissions) for the multimaterial beams. Minimizes compliance subject to a volume constraint
-  `mm.jl`: minimizes $CO_2$ emissions subject to compliance and volume fraction constraints
-  `mm_alternative.jl`: minimizes $CO_2$ emissions subject to displacement and mass constraints
- `post_fea_sm.jl` and `post_fea_mm.jl`: reruns the FEA of the final optimized designs, saving more detailed strain fields

The single- and multimaterial problems are solved for other geometries in `bridge/` and `Lbeam/`.

To run the examples, start Julia from the root directory with:
```console
// sync environment from Manifest.toml
$ julia --project=. -e "using Pkg; Pkg.instantiate()"

$ julia --project=. -t auto
julia> include("mm.jl")
julia> mm(100, 0.5, 5.0, 7.0, filename="path/to/results")
```

> [!TIP]
> If `filename` is set to `"results/mbb"`, the script will create the files:
> - `results/mbb.xxxx.vtu`: VTK Unstructured Grid File for each iteration
> - `results/mbb.pvd`: Paraview Data File to visualize all iterations as a time series
> - `results/mbb.h5`: HDF5 file containing the convergence history. Data is flattened to be compatible with other HDF5 readers
> - `results/mbb.xxxx.jld2`: JLD2 file containing the design variables at the final iteration. It is HDF5-compatible but best suited for post-processing by Julia code
>
> The full path to the parent directory (`results/` in this case) must already exist, and existing files in it may be overwritten
