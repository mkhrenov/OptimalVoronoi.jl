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