function lloyd(points::DMT, Ω::F; SMT=SparseMatrixCSC{Int,Int}, tol=1e-3, max_iters=1000) where {DMT,F}
    n = size(points, 2)
    z₀ = vec(points)

    f(x) = eval_cvt_objective(x, Ω)
    g!(g, x) = eval_cvt_objective_gradient!(g, x, Ω)

    o = BFGSOptimizerState(f, g!, Ω, length(z₀); f_abs_tol=1e-8)

    optimize_voronoi!(o, z₀)

    optimal_points = reshape(o.xₖ₊₁, 3, :)

    return bounded_voronoi(optimal_points, Ω)
end

function eval_cvt_objective(x, Ω::F) where {F}
    points = reshape(x, 3, :)
    voronoi = bounded_voronoi(points, Ω)

    cumulative_squared_distance(voronoi) / 10.0
end

function eval_cvt_objective_gradient!(g, x, Ω::F) where {F}
    points = reshape(x, 3, :)
    grad_p = reshape(g, 3, :)
    voronoi = bounded_voronoi(points, Ω)

    cumulative_squared_distance_gradient!(grad_p, voronoi)

    grad_p ./= 10.0
end

function cumulative_squared_distance(complex::CellComplex)
    J = 0.0
    for i in 1:n_cells(complex)
        z_i = SVector{3}(complex.cell_centers[1, i], complex.cell_centers[2, i], complex.cell_centers[3, i])
        f(y) = dot(y - z_i, y - z_i)

        J += cell_volume_integral(complex, f, i)
    end

    return J
end

function cumulative_squared_distance_gradient!(grad_p, complex::CellComplex)
    complex_centroids!(grad_p, complex)
    grad_p .*= -1.0
    grad_p .+= complex.cell_centers

    for i in 1:n_cells(complex)
        grad_p[:, i] .*= 2.0 * cell_volume(complex, i)
    end
end