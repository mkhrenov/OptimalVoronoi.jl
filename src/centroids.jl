

"""
    Calculate centroids of Voronoi cells used density function `f(x)`
"""
function get_centroids!(domain::AbstractArray, points::AbstractMatrix, N::Int, f::F) where {F}
    dims = size(domain)
    D = length(dims)
    # points = zeros(D, N)
    points .= 0
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

function get_centroids!(domain::AbstractArray, points::AbstractMatrix, N::Int)
    f(x) = 1.0
    get_centroids!(domain, points, N, f)
end
