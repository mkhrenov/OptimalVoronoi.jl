using OptimalVoronoi
using SparseArrays
using GLMakie
using Profile
using Random

sdf = ones(30, 40, 30) * Inf
# OptimalVoronoi.sdf_box!(sdf, 15, 20, 10, 24, 36, 10)
sdf_sphere!(sdf, 15, 16, 15, 8)
sdf_sphere!(sdf, 15, 18, 15, 8)
sdf_sphere!(sdf, 15, 20, 15, 8)
sdf_sphere!(sdf, 15, 22, 15, 8)
sdf_sphere!(sdf, 15, 24, 15, 8)

Ω = Ω_from_array(sdf)

points = sample_from_discrete_sdf(sdf, 100)

pcopy = copy(points)
scatter(points)

points .= pcopy
@time voronoi = centroidal_voronoi(points, Ω);
voronoi.vertices .-= 0.5
voronoi.cell_centers .-= 0.5

ρ = 7000.0e-9
k = 30.0e-3
cₚ = 500.0
h = 10.0e-6

A, e = mesh_fv_matrix_vector(voronoi, ρ, k, cₚ, h)

dt = 0.01
Ad = I + dt * A
ed = dt * e

vals = zeros(100)
vals[1] = 1.0

for k in 1:100
    vals .= Ad * vals + ed * 0.0
end

viz(voronoi, cell_colors=vals)
# volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)