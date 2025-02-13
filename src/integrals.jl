
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

    for f in faces_in_cell
        vol += frustum_volume(complex, cell, f)
    end

    return vol
end

function cell_centroid(complex::CellComplex{DMT,SMT}, cell::Int) where {DMT,SMT}
    faces_in_cell = view(rowvals(complex.E2T), nzrange(complex.E2T, cell))
    cell_centroid = zero(SVector{3})
    vol = 0.0

    for f in faces_in_cell
        frust_centroid, frust_vol = frustum_centroid(complex, cell, f)
        cell_centroid += frust_centroid * frust_vol
        vol += frust_vol
    end

    return cell_centroid / vol, vol
end

function frustum_volume(complex::CellComplex{DMT,SMT}, cell::Int, face::Int) where {DMT,SMT}
    # fc = project_to_face(complex, face, @view complex.cell_centers[:, cell])
    fc = face_centroid(complex, face)
    vol = 0.0

    edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, f))

    for e in edges_in_face
        edge_vertices = collect(complex.edge_vertex_set[e])

        vol += simplex_volume(
            view(complex.vertices, :, edge_vertices[1]),
            view(complex.vertices, :, edge_vertices[2]),
            fc,
            view(complex.cell_centers, :, cell)
        )
    end

    return vol
end

function frustum_centroid(complex::CellComplex{DMT,SMT}, cell::Int, face::Int) where {DMT,SMT}
    frust_centroid = zero(SVector{3})
    # fc = project_to_face(complex, face, @view complex.cell_centers[:, cell])
    fc = face_centroid(complex, face)
    vol = 0.0

    edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, face))

    for e in edges_in_face
        edge_vertices = collect(complex.edge_vertex_set[e])

        simp_vol = simplex_volume(
            view(complex.vertices, :, edge_vertices[1]),
            view(complex.vertices, :, edge_vertices[2]),
            fc,
            view(complex.cell_centers, :, cell)
        )
        frust_centroid += simplex_centroid(
            view(complex.vertices, :, edge_vertices[1]),
            view(complex.vertices, :, edge_vertices[2]),
            fc,
            view(complex.cell_centers, :, cell)
        ) * simp_vol
        vol += simp_vol
    end

    return frust_centroid / vol, vol
end

function simplex_volume(a, b, c, d)
    av = SVector{3}(a[1], a[2], a[3])
    bv = SVector{3}(b[1], b[2], b[3])
    cv = SVector{3}(c[1], c[2], c[3])
    dv = SVector{3}(d[1], d[2], d[3])

    return (1 / 6) * abs(det(
        hcat(
            av - dv, bv - dv, cv - dv
        )
    ))
end

function simplex_centroid(a, b, c, d)
    av = SVector{3}(a[1], a[2], a[3])
    bv = SVector{3}(b[1], b[2], b[3])
    cv = SVector{3}(c[1], c[2], c[3])
    dv = SVector{3}(d[1], d[2], d[3])

    return (av + bv + cv + dv) / 4.0
end

function face_centroid(complex::CellComplex{DMT,SMT}, face) where {DMT,SMT}
    verts = complex.vertices
    face_centr = zero(SVector{3})
    area = 0.0

    triangle = zeros(3, 3)

    edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, face))

    v0 = view(rowvals(complex.E0T), nzrange(complex.E0T, edges_in_face[1]))[1]
    p0 = SVector{3}(verts[1, v0], verts[2, v0], verts[3, v0])

    triangle[:, 1] .= @view verts[:, v0]


    for e in @view edges_in_face[2:end]
        v = view(rowvals(complex.E0T), nzrange(complex.E0T, e))
        p1 = SVector{3}(verts[1, v[1]], verts[2, v[1]], verts[3, v[1]])
        p2 = SVector{3}(verts[1, v[2]], verts[2, v[2]], verts[3, v[2]])

        triangle_centroid = (p0 + p1 + p2) / 3.0
        triangle_area = (1 / 2) * norm((p1 - p0) × (p2 - p0))

        face_centr += triangle_centroid * triangle_area
        area += triangle_area
    end

    return face_centr / area
end

function project_to_face(complex::CellComplex{DMT,SMT}, face::Int, point) where {DMT,SMT}
    verts = complex.vertices
    fv = collect(complex.face_vertex_set[face])

    p0 = SVector{3}(verts[1, fv[1]], verts[2, fv[1]], verts[3, fv[1]])
    p1 = SVector{3}(verts[1, fv[2]], verts[2, fv[2]], verts[3, fv[2]])
    p2 = SVector{3}(verts[1, fv[3]], verts[2, fv[3]], verts[3, fv[3]])

    normal = (p1 - p0) × (p2 - p0)
    normal /= norm(normal)

    return point - (normal ⋅ point) * normal
end

# Integrals of arbitrary functions

