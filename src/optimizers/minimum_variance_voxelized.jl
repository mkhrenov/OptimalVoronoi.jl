function minimum_variance_voronoi_vox(points::DMT, Ω::F, u::G, domain, sqdist, volumes, A, e; tol=1e-3, max_iters=1000) where {DMT,F,G}
    n = size(points, 2)
    z₀ = vec(points)
    ū = cu(zeros(1, n))

    f(z) = eval_minvar_vox_objective(z, domain, sqdist, u, ū, volumes)
    g!(g, z) = eval_minvar_vox_gradient!(g, z, domain, A, e, u, ū, volumes)

    o = BFGSOptimizerState(f, g!, Ω, length(z₀), size(points, 2); f_abs_tol=1e-6, f_rel_tol=1e-7, μ₀=1e-2, ST=Float32, VT=CuVector{Float32}, MT=CuMatrix{Float32}, max_iters=max_iters)

    optimize_voronoi!(o, z₀)

    optimal_points = reshape(o.xₖ₊₁, 3, :)

    return optimal_points
end

function eval_minvar_vox_objective(x, domain, sqdist, u::G, ū, volumes) where {G}
    points = reshape(x, 3, :)
    color_voronoi!(domain, points)
    cell_averages!(ū, volumes, domain, u)

    return (cumulative_squared_distance_vox(points, domain, sqdist) * 5e-6 + cumulative_variance_vox(domain, sqdist, u, ū)) * 1e-1
end

function eval_minvar_vox_gradient!(g, x, domain, A, e, u::G, ū, volumes) where {G}
    gm = reshape(g, 3, :)
    points = reshape(x, 3, :)
    color_voronoi!(domain, points)
    cell_averages!(ū, volumes, domain, u)

    cumulative_squared_distance_vox_gradient!(gm, points, domain)
    gm .*= 5e-6
    variance_vox_gradient!(gm, domain, points, u, ū, A, e)
    gm .*= 1e-1
end

# function cell_averages(complex::CellComplex{DMT,SMT}, u::G) where {DMT,SMT,G}
#     ū = DMT(zeros(1, n_cells(complex)))
#     cell_averages!(ū, complex, u)
#     return ū
# end

function cell_averages!(ū, volumes, domain, u::G) where {G}
    ū .= 0
    volumes .= 0
    f(p, i) = 1
    g(p, i) = u(p)

    cell_volume_integrals!(volumes, f, domain)
    cell_volume_integrals!(ū, g, domain)
    ū ./= volumes
end

function cumulative_variance_vox(domain, sqerr, u::G, ū) where {G}
    g(p, i) = (u(p) - ū[i])^2

    sqerr .= 0
    cell_volume_integrals!(sqerr, g, domain)

    return sum(sqerr)
end

function variance_vox_gradient!(gm, domain, points, u::G, ū, A, e) where {G}
    gx(p, i, j) = (ū[i]^2 - ū[j]^2 + 2(ū[j] - ū[i]) * u(p)) * (p[1] - points[1, i]) / sqrt((points[1, i] - points[1, j])^2 + (points[2, i] - points[2, j])^2 + (points[3, i] - points[3, j])^2)
    gy(p, i, j) = (ū[i]^2 - ū[j]^2 + 2(ū[j] - ū[i]) * u(p)) * (p[2] - points[2, i]) / sqrt((points[1, i] - points[1, j])^2 + (points[2, i] - points[2, j])^2 + (points[3, i] - points[3, j])^2)
    gz(p, i, j) = (ū[i]^2 - ū[j]^2 + 2(ū[j] - ū[i]) * u(p)) * (p[3] - points[3, i]) / sqrt((points[1, i] - points[1, j])^2 + (points[2, i] - points[2, j])^2 + (points[3, i] - points[3, j])^2)

    adjacency_matrix_vector!(A, e, domain)
    areas = sparse(Float32.(A))
    # o = CUDA.ones(size(areas, 2))

    nonzeros(areas) .= 0.0
    neighbor_surface_integrals!(areas, gx, domain, points)
    view(gm, 1, :) .+= vec(sum(areas, dims=1))

    nonzeros(areas) .= 0.0
    neighbor_surface_integrals!(areas, gy, domain, points)
    view(gm, 2, :) .+= vec(sum(areas, dims=1))

    nonzeros(areas) .= 0.0
    neighbor_surface_integrals!(areas, gz, domain, points)
    view(gm, 3, :) .+= vec(sum(areas, dims=1))
end
