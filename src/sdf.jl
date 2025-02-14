
# Treat indices as XYZ
function trilinear_interpolation(coord, array)
    Nx, Ny, Nz = size(array)
    x0 = floor(Int, coord[1])
    y0 = floor(Int, coord[2])
    z0 = floor(Int, coord[3])

    x1 = ceil(Int, coord[1])
    y1 = ceil(Int, coord[2])
    z1 = ceil(Int, coord[3])

    if (x1 > Nx || y1 > Ny || z1 > Nz ||
        x0 < 1 || y0 < 1 || z0 < 1)
        return 1.0
    end

    xd = (coord[1] - x0) / (x1 - x0)
    yd = (coord[2] - y0) / (y1 - y0)
    zd = (coord[3] - z0) / (z1 - z0)

    # Interpolate along X
    c00 = array[x0, y0, z0] * (1 - xd) + array[x1, y0, z0] * xd
    c01 = array[x0, y0, z1] * (1 - xd) + array[x1, y0, z1] * xd
    c10 = array[x0, y1, z0] * (1 - xd) + array[x1, y1, z0] * xd
    c11 = array[x0, y1, z1] * (1 - xd) + array[x1, y1, z1] * xd

    # Interpolate along Y
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    # Interpolate along Z
    return c0 * (1 - zd) + c1 * zd
end

function Ω_from_array(array)
    function Ω(x)
        return trilinear_interpolation(x, array)
    end
end

# 3 x N points
function in_Ω(points, Ω::F) where {F}
    return vec(mapslices(x -> Ω(x) < 0.0, points, dims=1))
end

function intersect_sdf(init, dir, Ω::F; tol=1e-2) where {F}
    t = 1.0
    # Use dir vector and start vertex to get intersection with ∂Ω via SDF
    for k in 1:20
        ω = Ω(init + t * dir)
        if abs(ω) < tol
            break
        end
        t -= ω
    end

    return init + t * dir
end

function sample_from_discrete_sdf(sdf, N)
    M = ndims(sdf)
    cidxs = CartesianIndices(sdf)
    points = [Float64(cidx[i]) for i in 1:M, cidx in shuffle(cidxs[(sdf).<(-0.5)])[1:N]]
    points .+= rand(size(points)...)
    points .-= 0.5

    return points
end