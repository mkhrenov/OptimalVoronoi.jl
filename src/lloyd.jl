function lloyd(points::DMT, Ω::F; SMT=SparseMatrixCSC, tol=1e-2, max_iters=1000) where {DMT,F}
    centroids = zeros(size(points))
    voronoi = bounded_voronoi(points, Ω; SMT=SMT)

    for k in 1:max_iters
        complex_centroids!(centroids, voronoi)

        if maximum(abs, points .- centroids) < tol
            break
        end

        points .= centroids
        voronoi = bounded_voronoi(points, Ω; SMT=SMT)
    end

    return voronoi
end