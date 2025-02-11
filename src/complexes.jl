# Incidence matrix representation for a 3-complex
struct CellComplex{DMT,SMT}
    vertices::DMT
    cell_centers::DMT

    E0::SMT # Vertex to edge incidence matrix
    E1::SMT # Edge to face incidence matrix
    E2::SMT # Face to volume incidence matrix

    E0T::SMT # Edge to vertex incidence matrix
    E1T::SMT # Face to edge incidence matrix
    E2T::SMT # Volume to face incidence matrix

    function CellComplex{DMT,SMT}(vertices::DMT, cell_centers::DMT, E0::SMT, E1::SMT, E2::SMT) where {DMT,SMT}
        E0T = SMT(transpose(E0))
        E1T = SMT(transpose(E1))
        E2T = SMT(transpose(E2))

        return new{DMT,SMT}(
            vertices, cell_centers,
            E0, E1, E2,
            E0T, E1T, E2T)
    end
end

n_vertices(complex::CellComplex) = size(complex.vertices, 1)
n_edges(complex::CellComplex) = size(complex.E0, 1)
n_faces(complex::CellComplex) = size(complex.E1, 1)
n_volumes(complex::CellComplex) = size(complex.E2, 1)

struct SubComplex{DMT,SMT,V}
    parent::CellComplex{DMT,SMT}

    vertex_sub::V
    edge_sub::V
    face_sub::V
    volume_sub::V
end


# Get a vector with 1s for each face that is on the boundary (is only incident to one volume)
function boundary_faces(complex::CellComplex)
    boundary = vec(sum(complex.E2, dims=1))

    return boundary .== 1
end

function get_boundary(complex::CellComplex)
    volume_sub = falses(n_volumes(complex))
    face_sub = boundary_faces(complex)
    edge_sub = (complex.E1T * face_sub) .> 0
    vertex_sub = (complex.E0T * edge_sub) .> 0

    return SubComplex(complex, vertex_sub, edge_sub, face_sub, volume_sub)
end

function prune_points(complex::CellComplex, vertices_to_drop)
    edges_to_drop = (complex.E0 * vertices_to_drop) .> 0
    faces_to_drop = (complex.E1 * edges_to_drop) .> 0
    volume_sub = trues(n_volumes(complex))

    return SubComplex(complex, .!vertices_to_drop, .!edges_to_drop, .!faces_to_drop, volume_sub)
end

function triangle_normal(p0, p1, p2, o)
    v1 = SVector{3}(p1[1], p1[2], p1[3]) - SVector{3}(p0[1], p0[2], p0[3])
    v2 = SVector{3}(p2[1], p2[2], p2[3]) - SVector{3}(p0[1], p0[2], p0[3])
    ov = SVector{3}(o[1], o[2], o[3])

    n = cross(v1, v2)

    if (ov - p0) ⋅ n > 0.0
        n = n * -1
    end

    return n / norm(n)
end