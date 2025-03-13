"""
    cell_to_voxel_map(domain)

    Takes a colored domain and returns a sparse matrix mapping each cell to voxel
"""
function cell_to_voxel_map(domain::Array{T,D}) where {T,D}
    cols = vec(domain) .+ 1
    rows = 1:length(domain)
    vals = ones(T, length(domain))

    return sparse(rows, cols, vals)
end

function cell_to_voxel_map(domain::CuArray{T,D,M}) where {T,D,M}
    cols = vec(domain) .+ 1
    rows = CuArray{T,1,M}(1:length(domain))
    vals = CuArray{T,1,M}(CUDA.ones(length(domain)))

    return sparse(rows, cols, vals)
end

function cell_volumes(domain, n_cells)
    volumes = zeros(n_cells)

    for i in domain
        if i == 0
            continue
        end

        volumes[i] += 1.0
    end

    return volumes
end

function cell_volumes(domain::CuArray{T,D,M}, n_cells) where {T,D,M}
    volumes = CUDA.zeros(n_cells)
    nthreads = 256
    nblocks = ceil(Int, length(domain) / nthreads)

    @cuda threads = nthreads blocks = nblocks volume_kernel!(domain, volumes)
    return volumes
end

function volume_kernel!(domain::CuDeviceArray{T1,D,M}, volumes::CuDeviceVector{T2,M}) where {T1,T2,D,M}
    i = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if i > length(domain)
        return nothing
    end

    if domain[i] != 0
        CUDA.@atomic volumes[domain[i]] += 1.0
    end

    return nothing
end

function cell_to_cell_map(domain1, domain2, n_cells_1, n_cells_2)
    incidence = zeros(n_cells_2, n_cells_1)

    for (i, j) in zip(domain2, domain1)
        if i == 0 || j == 0
            continue
        end

        incidence[i, j] += 1.0
    end

    incidence ./= cell_volumes(domain2, n_cells_2)
    return sparse(incidence)
end

function cell_to_cell_map(domain1::CuArray{T,D,M}, domain2::CuArray{T,D,M}, n_cells_1, n_cells_2) where {T,D,M}
    incidence = CUDA.zeros(n_cells_2, n_cells_1)

    nthreads = 256
    nblocks = ceil(Int, length(domain1) / nthreads)

    @cuda threads = nthreads blocks = nblocks incidence_kernel!(domain1, domain2, incidence)

    incidence ./= cell_volumes(domain2, n_cells_2)
    return sparse(incidence)
end

function incidence_kernel!(domain1::CuDeviceArray{T1,D,M}, domain2::CuDeviceArray{T1,D,M}, incidence::CuDeviceMatrix{T2,M}) where {T1,T2,D,M}
    i = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if i > length(domain1)
        return nothing
    end

    if domain1[i] != 0 && domain2[i] != 0
        CUDA.@atomic incidence[domain2[i], domain1[i]] += 1.0
    end

    return nothing
end

function cell_to_cell_map(domain1, domain2)
    n_cells_1 = maximum(domain1)
    n_cells_2 = maximum(domain2)
    return cell_to_cell_map(domain1, domain2, n_cells_1, n_cells_2)
end