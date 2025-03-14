# Primitives to build up bead geometry
function sdf_box!(sdf, px, py, pz, lx, ly, lz)
    map!((cidx, d_o) -> begin
            x, y, z = cidx[1], cidx[2], cidx[3]
            dx = abs(x - px) / (lx / 2)
            dy = abs(y - py) / (ly / 2)
            dz = abs(z - pz) / (lz / 2)

            d = 0.0
            if dx > dy && dx > dz
                d = (dx - 1.0) * (lx / 2)
            elseif dy > dz
                d = (dy - 1.0) * (ly / 2)
            else
                d = (dz - 1.0) * (lz / 2)
            end

            return min(d_o, d)
        end, sdf, CartesianIndices(sdf), sdf)
end

function sdf_sphere!(sdf, px, py, pz, r)
    map!((cidx, d_o) -> begin
            x, y, z = cidx[1], cidx[2], cidx[3]
            return min(d_o, √((x - px)^2 + (y - py)^2 + (z - pz)^2) - r)
        end, sdf, CartesianIndices(sdf), sdf)
end

# Treat indices as XYZ
function trilinear_interpolation(coord, array)
    Nx, Ny, Nz = size(array)
    x0 = floor(Int, coord[1])
    y0 = floor(Int, coord[2])
    z0 = floor(Int, coord[3])

    x1 = floor(Int, coord[1] + 1)
    y1 = floor(Int, coord[2] + 1)
    z1 = floor(Int, coord[3] + 1)

    if (x1 > Nx || y1 > Ny || z1 > Nz)
        return max(x1 - Nx, y1 - Ny, z1 - Nz)
    end
    if (x0 < 1 || y0 < 1 || z0 < 1)
        return max(1 - x0, 1 - y0, 1 - z0)
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
    points = [Float64(cidx[i]) for i in 1:M, cidx in shuffle(cidxs[(sdf).<(-3)])[1:N]]
    points .+= rand(size(points)...) .* 0.8
    points .-= 0.4

    return points
end

function eval_sdf!(Ω::G, c, points) where {G}
    map!(p -> Ω(p), c, eachcol(points))
end

function eval_sdf!(Ω::G, c::CuVector{T,M}, points::CuMatrix{T,M}) where {G,T,M}
    nthreads = 256
    nblocks = ceil(Int, length(c) / nthreads)

    @cuda threads = nthreads blocks = nblocks eval_sdf_kernel(Ω, c, points)
end

function eval_sdf_kernel(Ω::G, c::CuDeviceVector{T,M}, points::CuDeviceMatrix{T,M}) where {G,T,M}
    x = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x

    if x > length(c)
        return nothing
    end

    p = SVector{3}(points[1, x], points[2, x], points[3, x])
    c[x] = Ω(p)

    return nothing
end

function eval_sdf_grad!(Ω::G, g_points, points) where {G}
    for (gv, pv) in zip(eachcol(g_points), eachcol(points))
        finite_difference_gradient!(gv, Ω, pv)
    end
end

function eval_sdf_grad!(Ω::G, g_points::CuMatrix{T,M}, points::CuMatrix{T,M}) where {G,T,M}
    nthreads = 256
    nblocks = ceil(Int, size(points, 2) / nthreads)

    @cuda threads = nthreads blocks = nblocks eval_sdf_grad_kernel(Ω, g_points, points)
end

function eval_sdf_grad_kernel(Ω::G, g_points::CuDeviceMatrix{T,M}, points::CuDeviceMatrix{T,M}) where {G,T,M}
    x = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x

    if x > size(points, 2)
        return nothing
    end

    p = SVector{3}(points[1, x], points[2, x], points[3, x])
    finite_difference_gradient!(view(g_points, :, x), Ω, p)

    return nothing
end

function finite_difference_gradient!(g, f::F, x) where {F}
    h = 0.25
    dx = SVector{3}(h, 0.0, 0.0)
    dy = SVector{3}(0.0, h, 0.0)
    dz = SVector{3}(0.0, 0.0, h)

    g[1] = (f(x + dx) - f(x - dx)) / 2h
    g[2] = (f(x + dy) - f(x - dy)) / 2h
    g[3] = (f(x + dz) - f(x - dz)) / 2h
end