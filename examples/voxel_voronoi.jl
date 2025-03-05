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

points = sample_from_discrete_sdf(sdf, 500)

domain = map(x -> x ≤ 0.0 ? 1 : 0, sdf)

# Move data to GPU
domain = cu(domain)
points = cu(points)

@time color_voronoi!(domain, points)

# Move data back to CPU
domain = Array(domain)
points = Array(points)

fig = volume(domain)
scatter!(points, label=nothing)
display(fig)