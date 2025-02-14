
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

# # Integrals of arbitrary functions
# # function complex_centroids!(centroids::DMT, complex::CellComplex{DMT,SMT}) where {DMT,SMT}
# #     _, N = size(centroids)
# #     @assert N == n_volumes(complex)

# #     for i in 1:N
# #         centr, vol = cell_centroid(complex, i)
# #         centroids[:, i] .= centr
# #     end
# # end

# function cell_integral(complex::CellComplex{DMT,SMT}, Ω::F, α::K, cell::Int) where {DMT,SMT,F,K}
#     faces_in_cell = view(rowvals(complex.E2T), nzrange(complex.E2T, cell))
#     total = zero(0)

#     for face in faces_in_cell
#         total += frustum_integral(complex, Ω, α, cell, face)
#     end

#     return total
# end

# function frustum_integral(complex::CellComplex{DMT,SMT}, Ω::F, α::K, cell::Int, face::Int) where {DMT,SMT,F,K}
#     fc = face_centroid(complex, face)
#     total = zero(0)

#     edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, f))

#     for e in edges_in_face
#         edge_vertices = view(rowvals(complex.E0T), nzrange(complex.E0T, e))

#         total += simplex_integral(
#             view(complex.vertices, :, edge_vertices[1]),
#             view(complex.vertices, :, edge_vertices[2]),
#             fc,
#             view(complex.cell_centers, :, cell),
#             Ω, α
#         )
#     end

#     return total
# end

# function tetrahedron_integral(a, b, c, d, Ω::F, α::K) where {F,K}
#     av = SVector{3}(a[1], a[2], a[3])
#     bv = SVector{3}(b[1], b[2], b[3])
#     cv = SVector{3}(c[1], c[2], c[3])
#     dv = SVector{3}(d[1], d[2], d[3])

# end

# function triangle_integral(a, b, c, Ω::F, α::K) where {F,K}
#     av = SVector{3}(a[1], a[2], a[3])
#     bv = SVector{3}(b[1], b[2], b[3])
#     cv = SVector{3}(c[1], c[2], c[3])

#     u = av - cv
#     v = bv - cv

# end