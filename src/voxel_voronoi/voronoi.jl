
"""
    Color the N-dimensional `domain` (representing a regular voxel grid)
    based on which point in points it is closest too.

    Points in domain with value 0 are treated as out of bounds and will remain 0.
"""
function voronoi!(domain::AbstractArray, points::AbstractMatrix)
    dims = size(domain)
    indices = CartesianIndices(domain)
    D = size(points, 1)
    N = size(points, 2)

    for index in indices
        min_dist = typemax(Float64)

        for p in 1:N
            # Don't bother evaluating if the cell is out of bounds
            if domain[index] == 0
                continue
            end

            old_dist = min_dist
            new_dist = 0.0

            # Compute squared distance to point p
            for d in 1:D
                new_dist += (index[d] - points[d, p])^2
            end

            # Update if closer distance discovered
            if new_dist < old_dist
                min_dist = new_dist
                domain[index] = p
            end
        end
    end
end


"""
    Color the N-dimensional `domain` (representing a regular voxel grid)
    based on which point in points it is closest too.

    Points in domain with value 0 are treated as out of bounds and will remain 0.
"""
function voronoi!(domain::CuArray{T1,D,M}, points::CuMatrix{T2,M}) where {T1,T2,D,M}
    nthreads = 256
    nblocks = ceil(Int, length(domain) / nthreads)
    cidx = CartesianIndices(domain)

    @cuda threads = nthreads blocks = nblocks voronoi_kernel!(domain, points, cidx)
end

function voronoi_kernel!(domain::CuDeviceArray{T1,D,M}, points::CuDeviceMatrix{T2,M}, cidx::CartesianIndices) where {T1,T2,D,M}
    i = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if i > length(domain)
        return nothing
    end
    N = size(points, 2)

    index = cidx[i]
    min_dist = typemax(Float64)

    for p in 1:N
        # Don't bother evaluating if the cell is out of bounds
        if domain[index] == 0
            continue
        end

        old_dist = min_dist
        new_dist = 0.0

        # Compute squared distance to point p
        for d in 1:D
            new_dist += (index[d] - points[d, p])^2
        end

        # Update if closer distance discovered
        if new_dist < old_dist
            min_dist = new_dist
            domain[index] = p
        end
    end

    return nothing
end
