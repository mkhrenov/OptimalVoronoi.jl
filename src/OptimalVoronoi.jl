module OptimalVoronoi

using LinearAlgebra
using Random
using StaticArrays
using SparseArrays

using DataStructures
using ForwardDiff

using CUDA
using GLMakie
using Printf

export delaunay_tet, condense_delaunay, is_delaunay
export dual_complex, bound_voronoi, bounded_voronoi
export color_voronoi!, adjacency_matrix_vector!, cell_volume_integrals!, neighbor_surface_integrals!, cell_averages!, paint!

export Ω_from_array, in_Ω, sample_from_discrete_sdf
export sdf_box!, sdf_sphere!

export cell_averages, complex_volumes, complex_volume

export cell_to_cell_map, mesh_fv_matrix_vector

export centroidal_voronoi, minimum_variance_voronoi
export centroidal_voronoi_vox, minimum_variance_voronoi_vox
export viz!, viz

include("sdf.jl")
include("complexes/complexes.jl")

include("delaunay/delaunay.jl")

include("cell_complex_voronoi/cell_complex_voronoi.jl")
include("voxel_voronoi/voxel_voronoi.jl")

include("optimizers/optimizers.jl")

include("finite_volume/finite_volume.jl")

include("viz.jl")


end # module OptimalVoronoi
