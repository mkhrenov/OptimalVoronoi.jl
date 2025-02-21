using WeightedCVT
using SparseArrays
using GLMakie
using Profile
using Random

sdf = ones(30, 40, 30) * Inf
# WeightedCVT.sdf_box!(sdf, 15, 20, 10, 24, 36, 10)
WeightedCVT.sdf_sphere!(sdf, 15, 16, 15, 8)
WeightedCVT.sdf_sphere!(sdf, 15, 18, 15, 8)
WeightedCVT.sdf_sphere!(sdf, 15, 20, 15, 8)
WeightedCVT.sdf_sphere!(sdf, 15, 22, 15, 8)
WeightedCVT.sdf_sphere!(sdf, 15, 24, 15, 8)

Ω = WeightedCVT.Ω_from_array(sdf)

points = WeightedCVT.sample_from_discrete_sdf(sdf, 100)

pcopy = copy(points)
scatter(points)

points .= pcopy
@time voronoi = WeightedCVT.lloyd(points, Ω);

points .= pcopy
@time voronoi = WeightedCVT.lloyd(points, Ω);

# points .= pcopy
# @profview voronoi = WeightedCVT.lloyd(points, Ω);

# points .= pcopy
# @profview_allocs voronoi = WeightedCVT.lloyd(points, Ω) sample_rate=0.001

voronoi.vertices .-= 0.5
voronoi.cell_centers .-= 0.5
WeightedCVT.viz(voronoi)
volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)