
function adjacency_matrix_vector!(A::AbstractMatrix{T1}, e::AbstractVector{T2}, domain::AbstractArray{T3,D}) where {T1,T2,T3,D}
    A .= 0
    e .= 0

    cidxs = CartesianIndices(domain)

    for cidx in cidxs
        i = domain[cidx]

        if i == 0
            continue
        end

        for d in 1:D
            for offs in [-1, 1]
                offset = CartesianIndex(ntuple(i -> i == d ? offs : 0, D))
                nidx = cidx + offset

                if nidx in cidxs
                    j = domain[nidx]

                    if j == 0
                        e[i] = 1
                    elseif i != j
                        A[i, j] = 1
                    end
                end
            end
        end
    end
end

function adjacency_matrix_vector!(A::CuMatrix{T1,M}, e::CuVector{T2,M}, domain::CuArray{T3,D,M}) where {T1,T2,T3,D,M}
    A .= 0
    e .= 0

    nthreads = 256
    nblocks = ceil(Int, length(domain) / nthreads)
    cidxs = CartesianIndices(domain)

    @cuda threads = nthreads blocks = nblocks adjacency_kernel(A, e, domain, cidxs)
end

function adjacency_kernel(A::CuDeviceMatrix{T1,M}, e::CuDeviceVector{T2,M}, domain::CuDeviceArray{T3,D,M}, cidxs::CartesianIndices) where {T1,T2,T3,D,M}
    lidx = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if lidx > length(domain)
        return nothing
    end

    cidx = cidxs[lidx]
    i = domain[cidx]
    if i == 0
        return nothing
    end

    for d in 1:D
        for offs in -1:2:1
            offset = CartesianIndex(ntuple(i -> i == d ? offs : 0, D))
            nidx = cidx + offset

            if nidx in cidxs
                j = domain[nidx]

                if j == 0
                    e[i] = 1
                elseif i != j
                    A[i, j] = 1
                end
            end
        end
    end

    return nothing
end