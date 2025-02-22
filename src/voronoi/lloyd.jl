function lloyd(points::DMT, Ω::F; SMT=SparseMatrixCSC{Int,Int}, tol=1e-3, max_iters=1000) where {DMT,F}
    n = size(points, 2)
    z₀ = vec(points)

    f = eval_cvt_objective
    g! = eval_cvt_objective_gradient!
    c! = eval_cvt_constraint!

    optimal_points = reshape(optimize_voronoi(f, g!, c!, Ω, z₀, 3n, n; tol=tol, max_iters=max_iters), 3, :)

    return bounded_voronoi(optimal_points, Ω)
end


function eval_cvt_objective(complex)
    J = 0.0
    for i in 1:n_cells(voronoi)
        z_i = SVector{3}(complex.cell_centers[1, i], complex.cell_centers[2, i], complex.cell_centers[3, i])
        f(y) = dot(y - z_i, y - z_i)

        J += cell_volume_integral(complex, f, i)
    end

    return J / 100.0
end

function eval_cvt_objective_gradient!(grad_f, z, voronoi)
    points = reshape(z, 3, :)
    grad_p = reshape(grad_f, 3, :)

    complex_centroids!(grad_p, voronoi)
    grad_p .*= -1.0
    grad_p .+= points

    for i in 1:n_cells(voronoi)
        grad_p[:, i] .*= 2.0 * cell_volume(voronoi, i)
    end

    grad_p ./= 100.0
end

function eval_cvt_constraint!(c, z, Ω::F) where {F}
    points = reshape(z, 3, :)

    for i in 1:size(points, 2)
        z_i = SVector{3}(points[1, i], points[2, i], points[3, i])

        # Compute SDF distance violation
        c[i] = Ω(z_i)
    end
end