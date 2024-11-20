
"""
    Color the N-dimensional `domain` (representing a regular voxel grid)
    based on which point in points it is closest too.

    Points in domain with value 0 are treated as out of bounds and will remain 0.
"""
function voronoi!(domain::AbstractArray, min_dist::AbstractArray, points::AbstractMatrix)
    dims = size(domain)
    indices = CartesianIndices(domain)
    D = size(points, 1)
    N = size(points, 2)


    min_dist .= typemax(Float64)

    for index in indices
        for p in 1:N
            # Don't bother evaluating if the cell is out of bounds
            if domain[index] == 0
                continue
            end

            old_dist = min_dist[index]
            new_dist = 0.0

            # Compute squared distance to point p
            for d in 1:D
                new_dist += (index[d] - points[d, p])^2
            end

            # Update if closer distance discovered
            if new_dist < old_dist
                min_dist[index] = new_dist
                domain[index] = p
            end
        end
    end
end
