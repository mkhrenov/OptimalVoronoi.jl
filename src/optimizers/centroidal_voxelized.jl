function centroidal_voronoi_vox(points::DMT, Ω::F, domain, sqdist; tol=1e-3, max_iters=1000) where {DMT,F}
    n = size(points, 2)
    z₀ = vec(points)

    f(x) = eval_cvt_vox_objective(x, domain, sqdist)
    g!(g, x) = eval_cvt_vox_objective_gradient!(g, x, domain)

    o = BFGSOptimizerState(f, g!, Ω, length(z₀); f_abs_tol=1e-8, ST=Float32, VT=CuVector{Float32}, MT=CuMatrix{Float32})

    optimize_voronoi!(o, z₀)

    optimal_points = reshape(o.xₖ₊₁, 3, :)

    return optimal_points
end

function eval_cvt_vox_objective(x, domain, sqdist)
    points = reshape(x, 3, :)
    color_voronoi!(domain, points)

    cumulative_squared_distance_vox(points, domain, sqdist) / 10.0
end

function eval_cvt_vox_objective_gradient!(g, x, domain)
    points = reshape(x, 3, :)
    grad_p = reshape(g, 3, :)
    color_voronoi!(domain, points)

    cumulative_squared_distance_vox_gradient!(grad_p, points, domain)
end

function cumulative_squared_distance_vox(points, domain, sqdist)
    g(p, i) = (p[1] - points[1, i])^2 + (p[2] - points[2, i])^2 + (p[3] - points[3, i])^2
    sqdist .= 0
    cell_volume_integrals!(sqdist, g, domain)

    return sum(sqdist) / 1e8
end

function cumulative_squared_distance_vox_gradient!(grad_p, points, domain)
    weighted_centroid_diff = grad_p
    g(p, i) = SVector{3}(points[1, i], points[2, i], points[3, i]) - p

    cell_volume_integrals!(weighted_centroid_diff, g, domain)
    grad_p ./= 1e8
end