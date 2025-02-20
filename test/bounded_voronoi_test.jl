using WeightedCVT
using SparseArrays
using GLMakie

sdf = ones(20, 20, 20) * Inf
WeightedCVT.sdf_sphere!(sdf, 10, 10, 10, 8)
Ω = WeightedCVT.Ω_from_array(sdf)

points = 8 .* rand(3, 25) .+ 6

@time voronoi = WeightedCVT.bounded_voronoi(points, Ω);
voronoi.vertices .-= 0.5
WeightedCVT.viz(voronoi)
volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.05)