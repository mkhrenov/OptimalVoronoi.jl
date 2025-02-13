function average!(domain::AbstractArray, T::AbstractArray, volumes::AbstractVector, means::AbstractVector)
    @assert size(domain) == size(T)
    @assert length(volumes) == length(means)

    volumes .= 0
    means .= 0

    for index in CartesianIndices(domain)
        # If no membership, skip
        cell = domain[index]
        if cell == 0
            continue
        end

        volumes[cell] += 1
        means[cell] += T[index]
    end

    means ./= volumes

    return volumes, means
end

function average!(domain::CuArray, T::CuArray, volumes::CuVector, means::CuVector)
    @assert size(domain) == size(T)
    @assert length(volumes) == length(means)

    volumes .= 0
    means .= 0

    nthreads = 256
    nblocks = ceil(Int, length(domain) / nthreads)

    @cuda threads = nthreads blocks = nblocks average_kernel!(domain, T, volumes, means)

    means ./= volumes

    return volumes, means
end


function average_kernel!(domain::CuDeviceArray, T::CuDeviceArray, volumes::CuDeviceVector, means::CuDeviceVector)
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if i > length(domain)
        return nothing
    end

    # If no membership, skip
    cell = domain[index]
    if cell == 0
        return nothing
    end

    CUDA.@atomic volumes[cell] += 1
    CUDA.@atomic means[cell] += T[index]

    return nothing
end

function find_min_variance(points)
    t = delaunay_tet(points)
    dense_delaunay = condense_delaunay(t, points, SMT)
    dense_voronoi = dual_complex(dense_delaunay)
    
    return find_min_variance(points, dense_voronoi, dense_delaunay)
end

function find_min_variance(points, dense_voronoi, dense_delaunay)
    graph_delaunay = 0

    for i in 1:10


        # Check Delaunay condition to validate Voronoi topology after update
        # Only if it is violated, re-run the incremental flip algorithm
        if !is_delaunay(dense_delaunay, points)
            t = delaunay_tet(points)
            dense_delaunay = condense_delaunay(t, points, SMT)
            dense_voronoi = dual_complex(dense_delaunay)
        end
    end

    return dense_voronoi, dense_delaunay
end