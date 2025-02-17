function lloyd(points::DMT, Ω::F; SMT=SparseMatrixCSC{Int,Int}, tol=2e-2, max_iters=1000) where {DMT,F}
    centroids = zeros(size(points))
    next_points = copy(centroids)
    sdf = zeros(1, size(points, 2))
    α = ones(1, size(points, 2))
    voronoi = bounded_voronoi(points, Ω; SMT=SMT)

    for k in 1:max_iters
        complex_centroids!(centroids, voronoi)

        if maximum(abs, points .- centroids) < tol
            println("Lloyd took $k iterations")
            break
        end

        α .= 1.0
        @. next_points = α * (centroids - points) + points
        map!(Ω, sdf, eachslice(next_points, dims=2))
        while maximum(sdf) ≥ -0.01
            map!((x, y) -> x ≥ -0.01 ? (0.8 * y) : 1.0, α, sdf, α)
            @. next_points = α * (centroids - points) + points
            map!(Ω, sdf, eachslice(next_points, dims=2))
        end

        points .= next_points
        voronoi = bounded_voronoi(points, Ω; SMT=SMT)
    end

    return voronoi
end