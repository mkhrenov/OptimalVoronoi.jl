function minimum_variance_voronoi(points::DMT, Ω::F, u::G; tol=1e-3, max_iters=1000, α=1.0) where {DMT,F,G}
    n = size(points, 2)
    z₀ = vec(points)
    ū = zeros(1, n)

    f(z, Ω) = cumulative_error(z, Ω, u, ū)
    g!(g, z, Ω) = variance_gradient!(g, ū, z, Ω, u)
    c! = eval_cvt_constraint!

    optimal_points = reshape(optimize_voronoi(f, g!, c!, Ω, z₀, 3n, n; tol=tol, max_iters=max_iters, α=α), 3, :)

    return bounded_voronoi(optimal_points, Ω)
end

function variance_gradient!(gradient, ū, z, Ω::F, u::G) where {F,G}
    points = reshape(z, 3, :)
    voronoi = bounded_voronoi(points, Ω)
    cell_averages!(ū, voronoi, u)
    gradient .= 0
    gm = reshape(gradient, 3, :)

    complex_centroids!(gm, voronoi)
    gm .*= -1.0
    gm .+= points

    for i in 1:n_cells(voronoi)
        gm[:, i] .*= 2.0 * cell_volume(voronoi, i)
    end

    gm ./= 10.0

    for i in 1:n_cells(voronoi)
        for face in faces_of_cell(voronoi, i)
            for j in cells_of_face(voronoi, face)
                if i == j
                    continue
                end

                x_i = voronoi.cell_centers[:, i]
                x_j = voronoi.cell_centers[:, j]

                x_i = SVector{3}(x_i[1], x_i[2], x_i[3])
                x_j = SVector{3}(x_j[1], x_j[2], x_j[3])

                g(p) = (ū[i]^2 - ū[j]^2 + 2(ū[j] - ū[i]) * u(p)) * (p - x_i) / norm(x_i - x_j)

                gm[:, i] .+= face_surface_integral(voronoi, g, face)
            end
        end
    end

    gm ./= complex_volume(voronoi)
    gm .*= 10.0
end

function cell_averages!(ū::DMT, complex::CellComplex{DMT,SMT}, u::G) where {DMT,SMT,G}
    ū .= 0

    for cell in 1:n_cells(complex)
        ū[cell] = cell_volume_integral(complex, u, cell) / cell_volume(complex, cell)
    end
end

function cell_averages(complex::CellComplex{DMT,SMT}, u::G) where {DMT,SMT,G}
    ū = DMT(zeros(1, n_cells(complex)))
    cell_averages!(ū, complex, u)
    return ū
end

function cumulative_error(complex::CellComplex{DMT,SMT}, u::G) where {DMT,SMT,G}
    ū = DMT(zeros(1, n_cells(complex)))
    return cumulative_error(complex, u, ū)
end

function cumulative_error(z, Ω::F, u::G, ū) where {F,G}
    acc = 0.0
    points = reshape(z, 3, :)
    complex = bounded_voronoi(points, Ω)
    cell_averages!(ū, complex, u)

    for cell in 1:n_cells(complex)
        g(x) = (u(x) - ū[cell])^2
        acc += cell_volume_integral(complex, g, cell)
    end

    acc /= complex_volume(complex)

    acc += eval_cvt_objective(z, Ω) / 10.0

    return acc * 10.0
end