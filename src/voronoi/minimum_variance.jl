function minimum_variance_voronoi(points::DMT, Ω::F, u::G; tol=1e-3, max_iters=1000) where {DMT,F,G}
    n = size(points, 2)
    z₀ = vec(points)
    ū = zeros(1, n)

    f(z) = eval_minvar_objective(z, Ω, u, ū)
    g!(g, z) = eval_minvar_gradient!(g, z, Ω, u, ū)

    o = BFGSOptimizerState(f, g!, Ω, length(z₀); f_abs_tol=1e-8, μ₀=5e-3)

    optimize_voronoi!(o, z₀)

    optimal_points = reshape(o.xₖ₊₁, 3, :)

    return bounded_voronoi(optimal_points, Ω)
end

function eval_minvar_objective(x, Ω::F, u::G, ū) where {F,G}
    points = reshape(x, 3, :)
    voronoi = bounded_voronoi(points, Ω)

    return (5e-5) * cumulative_squared_distance(voronoi) + cumulative_variance(voronoi, u, ū)
end

function eval_minvar_gradient!(g, x, Ω::F, u::G, ū) where {F,G}
    gm = reshape(g, 3, :)
    points = reshape(x, 3, :)
    voronoi = bounded_voronoi(points, Ω)

    cumulative_squared_distance_gradient!(gm, voronoi)
    gm .*= 5e-5

    variance_gradient!(gm, voronoi, u, ū)
end

function cell_averages(complex::CellComplex{DMT,SMT}, u::G) where {DMT,SMT,G}
    ū = DMT(zeros(1, n_cells(complex)))
    cell_averages!(ū, complex, u)
    return ū
end

function cell_averages!(ū, complex::CellComplex{DMT,SMT}, u::G) where {DMT,SMT,G}
    ū .= 0
    map!(
        cell -> cell_volume_integral(complex, u, cell) / cell_volume(complex, cell),
        ū, 1:n_cells(complex)
    )
end

function cumulative_variance(complex::CellComplex, u::G, ū) where {G}
    cell_averages!(ū, complex, u)

    return sum(
        cell -> cell_volume_integral(complex, x -> (u(x) - ū[cell])^2, cell),
        1:n_cells(complex)
    )
end

function variance_gradient!(gm, complex::CellComplex, u::G, ū) where {G}
    cell_averages!(ū, complex, u)

    for i in 1:n_cells(complex)
        for face in faces_of_cell(complex, i)
            for j in cells_of_face(complex, face)
                if i == j
                    continue
                end

                x_i = complex.cell_centers[:, i]
                x_j = complex.cell_centers[:, j]

                x_i = SVector{3}(x_i[1], x_i[2], x_i[3])
                x_j = SVector{3}(x_j[1], x_j[2], x_j[3])

                g(p) = (ū[i]^2 - ū[j]^2 + 2(ū[j] - ū[i]) * u(p)) * (p - x_i) / norm(x_i - x_j)

                gm[:, i] .+= face_surface_integral(complex, g, face)
            end
        end
    end
end
