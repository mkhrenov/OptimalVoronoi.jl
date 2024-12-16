

"""
    Calculate centroids of Voronoi cells using density function `f(x)`
"""
function get_centroids!(domain::AbstractArray, points::AbstractMatrix, f::F) where {F}
    dims = size(domain)
    D = length(dims)
    points .= 0
    N = size(points, 2)
    masses = zeros(N)
    cur_point = zeros(D)

    for index in CartesianIndices(domain)
        # If no membership, skip
        cell = domain[index]
        if cell == 0
            continue
        end

        for d in 1:D
            cur_point[d] = index[d]
        end

        f_val = f(cur_point)
        masses[cell] += f_val
        for d in 1:D
            points[d, cell] += f_val * index[d]
        end
    end
    points ./= transpose(masses)

    return points
end

function get_centroids!(domain::AbstractArray, points::AbstractMatrix)
    f(x) = 1.0
    get_centroids!(domain, points, f)
end


"""
    Calculate centroids of Voronoi cells with uniform density
"""
function get_centroids!(domain::CuArray, points::CuMatrix)
    nthreads = 256
    nblocks = ceil(Int, length(domain) / nthreads)
    cidx = CartesianIndices(domain)

    points .= 0
    masses = CUDA.zeros(size(points, 2))

    @cuda threads = nthreads blocks = nblocks get_centroids_kernel!(domain, points, masses, cidx)

    points ./= transpose(masses)

    return points
end

function get_centroids_kernel!(domain::CuDeviceArray, points::CuDeviceMatrix, masses::CuDeviceVector, cidx::CartesianIndices)
    i = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if i > length(domain)
        return nothing
    end
    D = size(points, 1)

    index = cidx[i]

    # If no membership, skip
    cell = domain[index]
    if cell == 0
        return nothing
    end

    CUDA.@atomic masses[cell] += 1
    for d in 1:D
        CUDA.@atomic points[d, cell] += index[d]
    end

    return nothing
end