function cell_volume_integrals!(integrals::AbstractMatrix{T1}, f::F, domain::AbstractArray{T2,3}) where {F,T1,T2}
    cidxs = CartesianIndices(domain)

    for cidx in cidxs
        i = domain[cidx]

        if i == 0
            continue
        end

        x, y, z = cidx[1], cidx[2], cidx[3]
        p = SVector{3}(x, y, z)
        integrals[:, i] .+= f(p, i)
    end
end

function cell_volume_integrals!(integrals::CuMatrix{T1,M}, f::F, domain::CuArray{T2,3,M}) where {F,T1,T2,M}
    cidxs = CartesianIndices(domain)
    nthreads = 256
    nblocks = ceil(Int, length(domain) / nthreads)

    @cuda threads = nthreads blocks = nblocks cell_volume_integral_kernel(integrals, f, domain, cidxs)
end

function cell_volume_integral_kernel(integrals::CuDeviceMatrix{T1,M}, f::F, domain::CuDeviceArray{T2,3,M}, cidxs::CartesianIndices) where {F,T1,T2,M}
    lidx = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if lidx > length(domain)
        return nothing
    end

    D = size(integrals, 1)

    cidx = cidxs[lidx]

    i = domain[cidx]
    if i == 0
        return nothing
    end

    x, y, z = cidx[1], cidx[2], cidx[3]
    p = SVector{3}(x, y, z)
    r = f(p, i)
    @assert length(r) == D

    for d in 1:D
        CUDA.@atomic integrals[d, i] += r[d]
    end

    return nothing
end


"""
    `integrals` must already encode adjacency in its sparsity pattern
"""
function neighbor_surface_integrals!(integrals::AbstractSparseMatrix{T1,T2}, f::F, domain::AbstractArray{T3,3}, points::AbstractMatrix{T1}) where {F,T1,T2,T3}
    Ncells, _ = size(integrals)
    integral_vals = nonzeros(integrals)
    rvs = rowvals(integrals)
    e_z = SVector{3,T1}(0.0, 0.0, 1.0)

    for i in 1:Ncells
        p_i = SVector{3,T1}(points[1, i], points[2, i], points[3, i])

        for idx in nzrange(integrals, i)
            j = rvs[idx]
            p_j = SVector{3,T1}(points[1, j], points[2, j], points[3, j])
            p_mid = (p_i + p_j) / 2
            n = (p_j - p_i) / norm(p_j - p_i)

            # First in-plane basis vector
            u_1 = n × e_z
            u_1 /= norm(u_1)

            # Second in-plane basis vector
            u_2 = n × u_1

            # Rotation matrix
            R = hcat(n, u_1, u_2)

            # Sweep outwards in polar coordinates until a circle is fully beyond the face
            # (Beyond can be either in a third cell, or beyond the domain of integration)
            dr = 0.1
            for r in dr:dr:100.0
                any_valid = false
                dθ = min((2π / 10.0) / r / dr, 2π / 10.0)

                for θ in 0.0:dθ:2π
                    p_original = SVector{3,T1}(0.0, r * cos(θ), r * sin(θ))
                    p_plane = R * p_original + p_mid

                    # Check if point is within the face (is within the domain and is part of the boundary between cells i and j)
                    cidx = CartesianIndex(round(Int, p_plane[1]), round(Int, p_plane[2]), round(Int, p_plane[3]))
                    if !checkbounds(Bool, domain, cidx) || domain[cidx] == 0 || (domain[cidx] != i && domain[cidx] != j)
                        continue
                    end

                    any_valid = true
                    integral_vals[idx] += f(p_plane, i, j) * r * dr * dθ
                end

                # If none of the points in this radius were valid, we're done by assumption of convexity
                if !any_valid
                    break
                end
            end
        end
    end
end

function neighbor_surface_integrals!(integrals::CUSPARSE.CuSparseMatrixCSC{T1,T2}, f::F, domain::CuArray{T3,3,M}, points::CuMatrix{T1,M}) where {F,T1,T2,T3,M}
    n_cells = size(points, 2)
    max_incidence = mapreduce((x, y) -> x - y, max, @view(integrals.colPtr[2:end]), @view(integrals.colPtr[1:(end-1)]))

    nthreads_x = 64
    nthreads_y = 4
    nblocks_x = ceil(Int, n_cells / nthreads_x)
    nblocks_y = ceil(Int, max_incidence / nthreads_y)

    @cuda threads = (nthreads_x, nthreads_y) blocks = (nblocks_x, nblocks_y) neighbor_surface_integral_kernel(integrals, f, domain, points)
end

function neighbor_surface_integral_kernel(integrals::CUSPARSE.CuSparseDeviceMatrixCSC{T1,T2}, f::F, domain::CuDeviceArray{T3,3,M}, points::CuDeviceMatrix{T1,M}) where {F,T1,T2,T3,M}
    x = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    y = (blockIdx().y - Int32(1)) * blockDim().y + threadIdx().y
    e_z = SVector{3,T1}(0.0, 0.0, 1.0)

    i = x
    if i > size(points, 2) || (integrals.colPtr[i+1] - integrals.colPtr[i]) < y
        return nothing
    end

    idx = (integrals.colPtr[i]:(integrals.colPtr[i+1]-1))[y]
    j = integrals.rowVal[idx]

    p_i = SVector{3,T1}(points[1, i], points[2, i], points[3, i])
    p_j = SVector{3,T1}(points[1, j], points[2, j], points[3, j])

    p_mid = (p_i + p_j) / 2
    n = (p_j - p_i) / norm(p_j - p_i)

    # First in-plane basis vector
    u_1 = n × e_z
    # @cuprintln(norm(u_1))
    u_1 /= norm(u_1)

    # Second in-plane basis vector
    u_2 = n × u_1

    # Rotation matrix
    R = hcat(n, u_1, u_2)

    # Sweep outwards in polar coordinates until a circle is fully beyond the face
    # (Beyond can be either in a third cell, or beyond the domain of integration)
    dr = 0.1
    for r in dr:dr:100.0
        any_valid = false
        dθ = min((2π / 10.0) / r / dr, 2π / 10.0)

        for θ in 0.0:dθ:2π
            p_original = SVector{3,T1}(0.0, r * cos(θ), r * sin(θ))
            p_plane = R * p_original + p_mid

            # Check if point is within the face (is within the domain and is part of the boundary between cells i and j)
            cidx = CartesianIndex(round(Int32, p_plane[1]), round(Int32, p_plane[2]), round(Int32, p_plane[3]))

            if !checkbounds(Bool, domain, cidx) || domain[cidx] == 0 || (domain[cidx] != i && domain[cidx] != j)
                continue
            end

            any_valid = true
            val = f(p_plane, i, j) * r * dr * dθ
            CUDA.@atomic integrals.nzVal[idx] += val
        end

        # If none of the points in this radius were valid, we're done by assumption of convexity
        if !any_valid
            break
        end
    end

    return nothing
end


"""
    Ray-cast
"""
function exterior_surface_integrals()

end