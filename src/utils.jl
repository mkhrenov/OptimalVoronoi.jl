using StaticArrays
using LinearAlgebra

function circumcenter(points::SMatrix{N, D, T, L}) where {N, D, T, L}
    C = zeros(MMatrix{D + 2, D + 2, T})

    for i in 1:(D+2)
        for j in 1:(D+2)
            if i == j
                C[i, j] = zero(T)
            elseif (i == 1) || (j == 1)
                C[i, j] = one(T)
            else
                C[i, j] = norm(points[i-1, :] - points[j-1, :])
            end
        end
    end

    M = -2 * inv(C)
    m = @view M[1, 2:end]
    c = transpose(points) * m / sum(m)

    return c
end
