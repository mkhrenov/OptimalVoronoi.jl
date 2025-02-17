
# Basic integrals (uniform density / scalar field)

function complex_centroids!(centroids::DMT, complex::CellComplex{DMT,SMT}) where {DMT,SMT}
    _, N = size(centroids)
    @assert N == n_cells(complex)

    for i in 1:N
        centr, vol = cell_centroid(complex, i)
        centroids[:, i] .= centr
    end
end

function cell_volume(complex::CellComplex{DMT,SMT}, cell::Int) where {DMT,SMT}
    vol = 0.0

    for face in faces_of_cell(complex, cell)
        vol += cone_volume(complex, cell, face)
    end

    return vol
end

function cell_centroid(complex::CellComplex{DMT,SMT}, cell::Int) where {DMT,SMT}
    cell_centroid = zero(SVector{3,Float64})
    vol = 0.0

    rv_f = rowvals(complex.E2T)
    for face in faces_of_cell(complex, cell)
        frust_centroid, frust_vol = cone_centroid(complex, cell, face)

        cell_centroid += frust_centroid * frust_vol
        vol += frust_vol
    end

    return cell_centroid / vol, vol
end

function cone_volume(complex::CellComplex{DMT,SMT}, cell::Int, face::Int) where {DMT,SMT}
    fc = face_centroid(complex, face)
    vol = 0.0

    for e in edges_of_face(complex, face)
        edge_vertices = verts_of_edge(complex, edge)
        e1::Int = edge_vertices[1]
        e2::Int = edge_vertices[2]

        vol += simplex_volume(
            view(complex.vertices, :, e1),
            view(complex.vertices, :, e2),
            fc,
            view(complex.cell_centers, :, cell)
        )
    end

    return vol
end

function cone_centroid(complex::CellComplex{DMT,SMT}, cell::Int, face::Int) where {DMT,SMT}
    con_centroid = zero(SVector{3,Float64})
    fc = face_centroid(complex, face)
    vol = 0.0

    for edge in edges_of_face(complex, face)
        edge_vertices = verts_of_edge(complex, edge)
        v1::Int = edge_vertices[1]
        v2::Int = edge_vertices[2]

        simp_vol = simplex_volume(
            view(complex.vertices, :, v1),
            view(complex.vertices, :, v2),
            fc,
            view(complex.cell_centers, :, cell)
        )
        con_centroid += simplex_centroid(
            view(complex.vertices, :, v1),
            view(complex.vertices, :, v2),
            fc,
            view(complex.cell_centers, :, cell)
        ) * simp_vol
        vol += simp_vol
    end

    return con_centroid / vol, vol
end

function simplex_volume(a, b, c, d)
    av = SVector{3,Float64}(a[1], a[2], a[3])
    bv = SVector{3,Float64}(b[1], b[2], b[3])
    cv = SVector{3,Float64}(c[1], c[2], c[3])
    dv = SVector{3,Float64}(d[1], d[2], d[3])

    return (1.0 / 6.0) * abs(det(
        hcat(
            av - dv, bv - dv, cv - dv
        )
    ))
end

function simplex_centroid(a, b, c, d)
    av = SVector{3,Float64}(a[1], a[2], a[3])
    bv = SVector{3,Float64}(b[1], b[2], b[3])
    cv = SVector{3,Float64}(c[1], c[2], c[3])
    dv = SVector{3,Float64}(d[1], d[2], d[3])

    return (av + bv + cv + dv) / 4.0
end

function face_centroid(complex::CellComplex{DMT,SMT}, face::Int) where {DMT,SMT}
    verts = complex.vertices
    face_centr = zero(SVector{3,Float64})
    area = 0.0

    edges_in_face = edges_of_face(complex, face)
    e1::Int = edges_in_face[1]

    v0::Int = verts_of_edge(complex, e1)[1]
    p0 = SVector{3,Float64}(verts[1, v0], verts[2, v0], verts[3, v0])

    for (i, edge) in enumerate(edges_in_face)
        if i == 1
            continue
        end
        vs = verts_of_edge(complex, edge)
        v1::Int = vs[1]
        v2::Int = vs[2]

        p1 = SVector{3,Float64}(verts[1, v1], verts[2, v1], verts[3, v1])
        p2 = SVector{3,Float64}(verts[1, v2], verts[2, v2], verts[3, v2])

        triangle_centroid = (p0 + p1 + p2) / 3.0
        triangle_area = (1 / 2) * norm((p1 - p0) × (p2 - p0))

        face_centr += triangle_centroid * triangle_area
        area += triangle_area
    end

    return face_centr / area
end

# Integrals of arbitrary functions

function cell_volume_integral(complex::CellComplex{DMT,SMT}, α::K, cell::Int) where {DMT,SMT,K}
    faces_in_cell = view(rowvals(complex.E2T), nzrange(complex.E2T, cell))
    total = zero(α(view(complex.cell_centers, :, cell)))

    for face in faces_in_cell
        total += cone_volume_integral(complex, α, cell, face)
    end

    return total
end

function cone_volume_integral(complex::CellComplex{DMT,SMT}, α::K, cell::Int, face::Int) where {DMT,SMT,K}
    fc = face_centroid(complex, face)
    total = zero(α(view(complex.cell_centers, :, cell)))

    edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, f))

    for e in edges_in_face
        edge_vertices = view(rowvals(complex.E0T), nzrange(complex.E0T, e))

        total += tetrahedron_volume_integral(
            view(complex.vertices, :, edge_vertices[1]),
            view(complex.vertices, :, edge_vertices[2]),
            fc,
            view(complex.cell_centers, :, cell),
            α
        )
    end

    return total
end

function tetrahedron_volume_integral(a, b, c, d, α::K; N=10) where {K}
    a_v = SVector{3}(a[1], a[2], a[3])
    b_v = SVector{3}(b[1], b[2], b[3])
    c_v = SVector{3}(c[1], c[2], c[3])
    d_v = SVector{3}(d[1], d[2], d[3])

    v1 = a_v - d_v
    v2 = b_v - d_v
    v3 = c_v - d_v

    du = 1 / N
    dv = 1 / N
    dw = 1 / N
    dA = abs(det(hcat(v1, v2, v3))) / (N - 1) / N / (N + 1)

    total = zero(α(a_v))

    for u in 0:(N-1)
        for v in 1:(N-u)
            for w in 1:(N-u-v)
                p = u * du * v1 + v * dv * v2 + w * dw * v3 + d_v
                total += α(p) * dA
            end
        end
    end

    return total
end

function cell_surface_integral(complex::CellComplex{DMT,SMT}, α::K, cell::Int) where {DMT,SMT,K}
    faces_in_cell = view(rowvals(complex.E2T), nzrange(complex.E2T, cell))
    total = zero(α(view(complex.cell_centers, :, cell)))

    for face in faces_in_cell
        total += cone_surface_integral(complex, α, cell, face)
    end

    return total
end

function cone_surface_integral(complex::CellComplex{DMT,SMT}, α::K, cell::Int, face::Int) where {DMT,SMT,K}
    fc = face_centroid(complex, face)
    total = zero(α(view(complex.cell_centers, :, cell)))

    edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, f))

    for e in edges_in_face
        edge_vertices = view(rowvals(complex.E0T), nzrange(complex.E0T, e))

        total += triangle_surface_integral(
            view(complex.vertices, :, edge_vertices[1]),
            view(complex.vertices, :, edge_vertices[2]),
            fc,
            α
        )
    end

    return total
end

function triangle_surface_integral(a, b, c, α::K; N=10) where {K}
    av = SVector{3}(a[1], a[2], a[3])
    bv = SVector{3}(b[1], b[2], b[3])
    cv = SVector{3}(c[1], c[2], c[3])

    v1 = av - cv
    v2 = bv - cv

    du = 1 / N
    dv = 1 / N
    dA = norm(v1 × v2) / N / (N + 1)

    total = zero(α(a_v))

    for u in 0:(N-1)
        for v in 1:(N-u)
            p = u * du * v1 + v * dv * v2 + cv
            total += α(p) * dA
        end
    end

    return total
end