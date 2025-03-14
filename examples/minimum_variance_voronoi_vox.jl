using OptimalVoronoi
using SparseArrays
using GLMakie
using Profile
using Random
using CUDA


sdf = ones(300, 400, 200) * Inf
sdf_box!(sdf, 150, 200, 60, 240, 360, 100)
sdf_sphere!(sdf, 150, 160, 110, 80)
sdf_sphere!(sdf, 150, 180, 110, 80)
sdf_sphere!(sdf, 150, 200, 110, 80)
sdf_sphere!(sdf, 150, 220, 110, 80)
sdf_sphere!(sdf, 150, 240, 110, 80)

T(x) = exp(-((x[1] - 150)^2 + (x[2] - 240)^2 + (x[3] - 200)^2) / 50^2) + 1.0
Tarr = map(T, CartesianIndices(sdf)) .* (sdf .≤ 0)

N_cells = 1000
points = sample_from_discrete_sdf(sdf, N_cells)

domain = map(x -> x ≤ 0.0 ? 1 : 0, sdf)

# Move data to GPU
sdf = cu(sdf)
Ω = Ω_from_array(sdf)
domain = cu(domain)
points = cu(points)
volumes = CUDA.zeros(1, N_cells)
centroids = copy(points)#CUDA.zeros(3, N_cells)
sqdist = CUDA.zeros(1, N_cells)
A = CUDA.zeros(Int, N_cells, N_cells)
e = CUDA.zeros(Int, N_cells)

###################### Needs constraints to be re-introduced, line-search re-enabled ################################
###################### Ideally should use make it so that so we don't need to compute the objective, prevent GPU synchronizations ###########
###################### Could use a penalty method rather than barrier / IP now that we're unshackled from cell complex construction again? ########

@time voronoi = OptimalVoronoi.centroidal_voronoi_vox(points, Ω, domain, sqdist);
points .= voronoi

@time voronoi = OptimalVoronoi.minimum_variance_voronoi_vox(points, Ω, T, domain, sqdist, volumes, A, e; max_iters=500);
points .= voronoi


# points .= pcopy
# @time voronoi = minimum_variance_voronoi(points, Ω, T);

# points .= pcopy
# @profview voronoi = minimum_variance_voronoi(points, Ω, T);

# points .= pcopy
# @profview_allocs voronoi = minimum_variance_voronoi(points, Ω, T) sample_rate=0.001

# voronoi.vertices .-= 0.5
# voronoi.cell_centers .-= 0.5

# viz(voronoi, cell_colors=vec(cell_averages(voronoi, T)))
# volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)



color_voronoi!(domain, points)

# Move data back to CPU
domain = Array(domain)
points = Array(points)
volumes = Array(volumes)
# centroids = Array(centroids)

colormap = collect(cgrad(:default))
colormap[1] = RGBAf(0, 0, 0, 0)
fig = volume(domain, colormap=colormap)
scatter!(points, label=nothing)
# scatter!(centroids, label=nothing)
display(fig)

Tarr2 = copy(Tarr)
cell_vals = copy(volumes)
OptimalVoronoi.cell_averages!(cell_vals, volumes, domain, T)
OptimalVoronoi.paint!(Tarr2, domain, cell_vals)

volume(Tarr, colormap=colormap, colorrange=(0.98,  maximum(Tarr)))
volume(Tarr2, colormap=colormap, colorrange=(0.98, maximum(Tarr)))
# scatter!(points, label=nothing)


# viz(voronoi, cell_colors=vec(complex_volumes(voronoi)))
# volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)
