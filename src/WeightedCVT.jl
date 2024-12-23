module WeightedCVT

using LinearAlgebra
using GeometryBasics
using Meshing
using StaticArrays
using SparseArrays
using CUDA

include("voxel_voronoi/voronoi.jl")
include("voxel_voronoi/centroids.jl")
include("voxel_voronoi/lloyd.jl")

include("structure.jl")
include("minimum_variance.jl")

end # module WeightedCVT
