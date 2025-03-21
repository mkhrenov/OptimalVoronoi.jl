using OptimalVoronoi
using SparseArrays
using GLMakie
using Profile
using Random
using LinearAlgebra
using CUDA


f(p, i) = 1
g(p, i) = p

sdf = ones(300, 400, 200) * Inf
sdf_box!(sdf, 150, 200, 60, 240, 360, 100)
sdf_sphere!(sdf, 150, 160, 110, 80)
sdf_sphere!(sdf, 150, 180, 110, 80)
sdf_sphere!(sdf, 150, 200, 110, 80)
sdf_sphere!(sdf, 150, 220, 110, 80)
sdf_sphere!(sdf, 150, 240, 110, 80)

N_cells = 200
points = sample_from_discrete_sdf(sdf, N_cells)

domain = map(x -> x ≤ 0.0 ? 1 : 0, sdf)

sdf = cu(sdf)
Ω = Ω_from_array(sdf)

# Move data to GPU
domain = cu(domain)
points = cu(points)
volumes = CUDA.zeros(1, N_cells)
centroids = copy(points)#CUDA.zeros(3, N_cells)
sqdist = CUDA.zeros(1, N_cells)

pcopy = copy(points)


points .= cu(pcopy)
@time voronoi = centroidal_voronoi_vox(points, Ω, domain, sqdist);

# points .= pcopy
# @time voronoi = centroidal_voronoi_vox(points, Ω, domain, sqdist);

# points .= cu(pcopy)
# @profview voronoi = centroidal_voronoi_vox(points, Ω, domain, sqdist);
points .= voronoi
# points .= pcopy
# @profview_allocs voronoi = centroidal_voronoi(points, Ω) sample_rate=0.001

# Adj, e = OptimalVoronoi.mesh_adjacency_matrix_vector(points, domain, 1e-3)
# Deg = Diagonal(vec(sum(Adj, dims=2)))

A, e = mesh_fv_matrix_vector(points, domain, Ω, 1e-3, 7000.0, 30.0, 500.0, 10.0)

color_voronoi!(domain, points)

# Move data back to CPU
domain = Array(domain)
points = Array(points)
volumes = Array(volumes)
# centroids = Array(centroids)

colormap = to_colormap(:viridis)
colormap[1] = RGBAf(0, 0, 0, 0)
fig = volume(domain, colormap=colormap)
scatter!(points, label=nothing)
# scatter!(centroids, label=nothing)
display(fig)

# voronoi.vertices .-= 0.5
# voronoi.cell_centers .-= 0.5
# viz(voronoi)
# volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)