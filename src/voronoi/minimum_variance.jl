function minimum_variance_voronoi(points::DMT, Ω::F, u::G; SMT=SparseMatrixCSC{Int,Int}, tol=1e-3, η=1.0, max_iters=1000) where {DMT,F,G}
    gradient = DMT(zeros(size(points)))
    centroids = DMT(zeros(size(points)))
    averages = DMT(zeros(1, size(points, 2)))
    next_points = copy(points)

    sdf = DMT(zeros(1, size(points, 2)))

    voronoi = bounded_voronoi(points, Ω; SMT=SMT)
    α = η

    err = cumulative_error(voronoi, u, averages)
    for k in 1:max_iters
        complex_centroids!(centroids, voronoi)
        variance_gradient!(gradient, averages, voronoi, u)
        gradient .+= (points .- centroids) .* 0.01
        @show err = cumulative_error(voronoi, u, averages)

        @show maximum(abs, gradient)
        if maximum(abs, gradient) < tol
            println("Minimum Variance Voronoi took $k iterations")
            break
        end

        @. next_points = points - α * gradient
        map!(Ω, sdf, eachslice(next_points, dims=2))
        while maximum(sdf) ≥ -0.01 #|| err - cumulative_error(bounded_voronoi(next_points, Ω; SMT=SMT), u) < dot(gradient, gradient) * α * 1e-5
            α *= 0.8
            @. next_points = points - α * gradient
            map!(Ω, sdf, eachslice(next_points, dims=2))
        end
        α = clamp(α * 1.2, 1e-12, η)

        @show maximum(sdf)
        @show α

        points .= next_points
        voronoi = bounded_voronoi(points, Ω; SMT=SMT)
        if α < 1e-11
            break
        end
    end

    @show gradient

    return voronoi
end

function variance_gradient!(gradient::DMT, ū::DMT, voronoi::CellComplex{DMT,SMT}, u::G) where {DMT,SMT,G}
    cell_averages!(ū, voronoi, u)
    gradient .= 0

    for i in 1:n_cells(voronoi)
        # @show i
        for face in faces_of_cell(voronoi, i)
            for j in cells_of_face(voronoi, face)
                if i == j
                    continue
                end
                # @show j

                x_i = voronoi.cell_centers[:, i]
                x_j = voronoi.cell_centers[:, j]

                x_i = SVector{3}(x_i[1], x_i[2], x_i[3])
                x_j = SVector{3}(x_j[1], x_j[2], x_j[3])

                g(p) = (ū[i]^2 - ū[j]^2 + 2(ū[j] - ū[i]) * u(p)) * (p - x_i) / norm(x_i - x_j)

                gradient[:, i] .+= face_surface_integral(voronoi, g, face)
            end
        end
    end
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

function cumulative_error(complex::CellComplex{DMT,SMT}, u::G, ū::DMT) where {DMT,SMT,G}
    acc = 0.0
    cell_averages!(ū, complex, u)

    for cell in 1:n_cells(complex)
        g(x) = (u(x) - ū[cell])^2
        acc += cell_volume_integral(complex, g, cell)
    end

    return acc
end