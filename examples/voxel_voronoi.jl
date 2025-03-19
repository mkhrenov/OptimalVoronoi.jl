using OptimalVoronoi
using GLMakie
using CUDA
using LinearAlgebra
using SparseArrays

sdf = ones(300, 400, 300) * Inf
# sdf_box!(sdf, 15, 20, 10, 24, 36, 10)
sdf_sphere!(sdf, 150, 160, 150, 80)
sdf_sphere!(sdf, 150, 180, 150, 80)
sdf_sphere!(sdf, 150, 200, 150, 80)
sdf_sphere!(sdf, 150, 220, 150, 80)
sdf_sphere!(sdf, 150, 240, 150, 80)


f(p, i) = 1
g(p, i) = p
h(p, i, j) = 1


N_cells = 500
points = sample_from_discrete_sdf(sdf, N_cells)

domain = map(x -> x ≤ 0.0 ? 1 : 0, sdf)

# Move data to GPU
domain = cu(domain)
points = cu(points)
A = CUDA.zeros(N_cells, N_cells)
e = CUDA.zeros(N_cells)

CUDA.@time color_voronoi!(domain, points)
CUDA.@time adjacency_matrix_vector!(A, e, domain)

volumes = CUDA.zeros(1, N_cells)
CUDA.@time cell_volume_integrals!(volumes, f, domain)

centroids = CUDA.zeros(3, N_cells)
CUDA.@time cell_volume_integrals!(centroids, g, domain)
centroids ./= volumes

areas = sparse(Float32.(A))
nonzeros(areas) .= 0.0
CUDA.@time neighbor_surface_integrals!(areas, h, domain, points)

# Move data back to CPU
domain = Array(domain)
points = Array(points)
A = Array(A)
e = Array(e)
volumes = Array(volumes)
centroids = Array(centroids)
areas = sparse(Array(areas))

colormap = to_colormap(:viridis)
colormap[1] = RGBAf(0, 0, 0, 0)
fig = volume(domain, colormap=colormap)
scatter!(points, label=nothing)
scatter!(centroids, label=nothing)
display(fig)

# domain_ext = map(x -> (x > 0) && (A[x, 6] > 0), domain)
# volume(domain_ext)