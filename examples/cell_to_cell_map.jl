using OptimalVoronoi
using GLMakie
using CUDA

sdf = ones(300, 400, 300) * Inf
# sdf_box!(sdf, 15, 20, 10, 24, 36, 10)
sdf_sphere!(sdf, 150, 160, 150, 80)
sdf_sphere!(sdf, 150, 180, 150, 80)
sdf_sphere!(sdf, 150, 200, 150, 80)
sdf_sphere!(sdf, 150, 220, 150, 80)
sdf_sphere!(sdf, 150, 240, 150, 80)

points1 = sample_from_discrete_sdf(sdf, 500)
points2 = sample_from_discrete_sdf(sdf, 500)

domain1 = map(x -> x ≤ 0.0 ? 1 : 0, sdf)
domain2 = map(x -> x ≤ 0.0 ? 1 : 0, sdf)

# Move data to GPU
domain1 = cu(domain1)
points1 = cu(points1)
domain2 = cu(domain2)
points2 = cu(points2)

@time color_voronoi!(domain1, points1)
@time color_voronoi!(domain2, points2)

@time c2c = OptimalVoronoi.cell_to_cell_map(domain1, domain2)

# Move data back to CPU
domain1 = Array(domain1)
points1 = Array(points1)
domain2 = Array(domain2)
points2 = Array(points2)

vals1 = cu(randn(500))

fig = volume(domain1)
scatter!(points1, label=nothing)
display(fig)

fig = volume(domain2)
scatter!(points2, label=nothing)
display(fig)

