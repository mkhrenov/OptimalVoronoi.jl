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
# WeightedCVT.sdf_box!(sdf, 15, 24, 15, 16, 16, 16)

Ω = WeightedCVT.Ω_from_array(sdf)
T(x) = (((x[1] - 15)^2 + (x[2] - 24)^2 + (x[3] - 24)^2) + 1e-4)^(-0.5) * exp(-((x[2] - 24) + (((x[1] - 15)^2 + (x[2] - 24)^2 + (x[3] - 24)^2))^(0.5)) / 5^2)

points = WeightedCVT.sample_from_discrete_sdf(sdf, 500)

pcopy = copy(points)

points .= pcopy
# voronoi = WeightedCVT.lloyd(points, Ω)
# points .= voronoi.cell_centers
@time voronoi = WeightedCVT.minimum_variance_voronoi(points, Ω, T; max_iters=1000);

# points .= pcopy
# @time voronoi = WeightedCVT.minimum_variance_voronoi(points, Ω, T);

# points .= pcopy
# @profview voronoi = WeightedCVT.minimum_variance_voronoi(points, Ω, T);

# points .= pcopy
# @profview_allocs voronoi = WeightedCVT.minimum_variance_voronoi(points, Ω, T) sample_rate=0.001
# @show WeightedCVT.cumulative_error(voronoi, T)

voronoi.vertices .-= 0.5
voronoi.cell_centers .-= 0.5


WeightedCVT.viz(voronoi, cell_colors=vec(WeightedCVT.cell_averages(voronoi, T)))
volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)

# WeightedCVT.viz(voronoi, cell_colors=vec(WeightedCVT.complex_volumes(voronoi)))
# volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)
