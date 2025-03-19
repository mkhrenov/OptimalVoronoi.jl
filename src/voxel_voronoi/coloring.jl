function paint!(dst, domain, cell_vals)
    map!(c -> c == 0 ? 0.0 : cell_vals[c], dst, domain)
end

"""
    Color the N-dimensional `domain` (representing a regular voxel grid)
    based on which point in points it is closest too.

    Points in domain with value 0 are treated as out of bounds and will remain 0.
"""
function color_voronoi!(domain::AbstractArray{T1,D}, points::AbstractMatrix{T2}) where {T1,T2,D}
    indices = CartesianIndices(domain)
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
function color_voronoi!(domain::CuArray{T1,D,M}, points::CuMatrix{T2,M}) where {T1,T2,D,M}
    nthreads = 256
    nblocks = ceil(Int, length(domain) / nthreads)
    cidx = CartesianIndices(domain)

    @cuda threads = nthreads blocks = nblocks voronoi_kernel!(domain, points, cidx)
end

function voronoi_kernel!(domain::CuDeviceArray{T1,D,M}, points::CuDeviceMatrix{T2,M}, cidx::CartesianIndices) where {T1,T2,D,M}
    i = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    # Don't bother evaluating if the cell is out of bounds
    if i > length(domain) || domain[i] == 0
        return nothing
    end

    N = size(points, 2)

    @inbounds index = cidx[i]
    min_dist = typemax(T2)

    for p in 1:N
        new_dist = zero(T2)

        # Compute squared distance to point p
        for d in 1:D
            @inbounds new_dist += (index[d] - points[d, p])^2
        end

        # Update if closer distance discovered
        if new_dist < min_dist
            min_dist = new_dist
            @inbounds domain[index] = p
        end
    end

    return nothing
end