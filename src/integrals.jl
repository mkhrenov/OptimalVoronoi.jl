
# Basic integrals (uniform density / scalar field)

function complex_centroids!(centroids::DMT, complex::CellComplex{DMT,SMT}) where {DMT,SMT}
    _, N = size(centroids)
    @assert N == n_volumes(complex)

    for i in 1:N
        centr, vol = cell_centroid(complex, i)
        centroids[:, i] .= centr
    end
end

function cell_volume(complex::CellComplex{DMT,SMT}, cell::Int) where {DMT,SMT}
    faces_in_cell = view(rowvals(complex.E2T), nzrange(complex.E2T, cell))
    vol = 0.0

    for face in faces_in_cell
        vol += frustum_volume(complex, cell, face)
    end

    return vol
end

function cell_centroid(complex::CellComplex{DMT,SMT}, cell::Int) where {DMT,SMT}
    cell_centroid = zero(SVector{3,Float64})
    vol = 0.0

    rv_f = rowvals(complex.E2T)
    for f_idx::Int in nzrange(complex.E2T, cell)
        face::Int = rv_f[f_idx]

        frust_centroid, frust_vol = frustum_centroid(complex, cell, face)

        cell_centroid += frust_centroid * frust_vol
        vol += frust_vol
    end

    return cell_centroid / vol, vol
end

function frustum_volume(complex::CellComplex{DMT,SMT}, cell::Int, face::Int) where {DMT,SMT}
    fc = face_centroid(complex, face)
    vol::Float64 = 0.0

    rv_e = rowvals(complex.E1T)
    rv_v = rowvals(complex.E0T)

    edges_in_face = nzrange(complex.E1T, face)
    for e_idx::Int in edges_in_face
        e::Int = rv_e[e_idx]

        edge_vertices = nzrange(complex.E0T, e)
        e1i::Int = edge_vertices[1]
        e2i::Int = edge_vertices[2]
        e1::Int = rv_v[e1i]
        e2::Int = rv_v[e2i]

        vol += simplex_volume(
            view(complex.vertices, :, e1),
            view(complex.vertices, :, e2),
            fc,
            view(complex.cell_centers, :, cell)
        )::Float64
    end

    return vol
end

function frustum_centroid(complex::CellComplex{DMT,SMT}, cell::Int, face::Int) where {DMT,SMT}
    frust_centroid::SVector{3,Float64} = zero(SVector{3,Float64})
    fc = face_centroid(complex, face)
    vol::Float64 = 0.0

    rv_e::Vector{Int} = rowvals(complex.E1T)
    rv_v::Vector{Int} = rowvals(complex.E0T)

    for e_idx in nzrange(complex.E1T, face)::UnitRange{Int}
        e::Int = rv_e[e_idx]
        edge_vertices::UnitRange{Int} = nzrange(complex.E0T, e)
        v1::Int = rv_v[edge_vertices[1]]
        v2::Int = rv_v[edge_vertices[2]]

        simp_vol = simplex_volume(
            view(complex.vertices, :, v1),
            view(complex.vertices, :, v2),
            fc,
            view(complex.cell_centers, :, cell)
        )
        frust_centroid += simplex_centroid(
            view(complex.vertices, :, v1),
            view(complex.vertices, :, v2),
            fc,
            view(complex.cell_centers, :, cell)
        ) * simp_vol
        vol += simp_vol
    end

    return frust_centroid / vol, vol
end

function simplex_volume(a, b, c, d)
    av = SVector{3,Float64}(a[1], a[2], a[3])
    bv = SVector{3,Float64}(b[1], b[2], b[3])
    cv = SVector{3,Float64}(c[1], c[2], c[3])
    dv = SVector{3,Float64}(d[1], d[2], d[3])

    return (1 / 6) * abs(det(
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
    face_centr::SVector{3,Float64} = zero(SVector{3,Float64})
    area::Float64 = 0.0

    rv_v::Vector{Int} = rowvals(complex.E0T)
    rv_e::Vector{Int} = rowvals(complex.E1T)

    edges_in_face::UnitRange{Int} = nzrange(complex.E1T, face)
    e1::Int = rv_e[edges_in_face[1]]

    v0i::Int = nzrange(complex.E0T, e1)[1]
    v0::Int = rv_v[v0i]
    p0 = SVector{3,Float64}(verts[1, v0], verts[2, v0], verts[3, v0])

    for (i, e_idx) in enumerate(edges_in_face)
        if i == 1
            continue
        end
        e::Int = rv_e[e_idx]
        v::UnitRange{Int} = nzrange(complex.E0T, e)
        v1::Int = rv_v[v[1]]
        v2::Int = rv_v[v[2]]
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

