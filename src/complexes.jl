# Incidence matrix representation for a 3-complex
struct CellComplex{DMT,SMT}
    vertices::DMT     # Geometric locations of vertices (3xN)
    cell_centers::DMT # Geometric locations of cell centers (3xN)

    E0::SMT # Vertex to edge incidence matrix
    E1::SMT # Edge to face incidence matrix
    E2::SMT # Face to volume incidence matrix

    E0T::SMT # Edge to vertex incidence matrix
    E1T::SMT # Face to edge incidence matrix
    E2T::SMT # Volume to face incidence matrix

    edge_vertex_set::Vector{Set{Int}}
    face_vertex_set::Vector{Set{Int}}
    cell_vertex_set::Vector{Set{Int}}

    function CellComplex{DMT,SMT}(vertices::DMT, cell_centers::DMT, E0::SMT, E1::SMT, E2::SMT) where {DMT,SMT}
        E0T = SMT(transpose(E0))
        E1T = SMT(transpose(E1))
        E2T = SMT(transpose(E2))

        N_edges = size(E0, 1)
        N_faces = size(E1, 1)
        N_cells = size(E2, 1)

        edge_vertex_set = Vector{Set{Int}}(undef, N_edges)
        face_vertex_set = Vector{Set{Int}}(undef, N_faces)
        cell_vertex_set = Vector{Set{Int}}(undef, N_cells)

        for e in 1:N_edges
            edge_vertex_set[e] = Set(rowvals(E0T)[nzrange(E0T, e)])
        end

        for f in 1:N_faces
            edges_in_face = view(rowvals(E1T), nzrange(E1T, f))
            face_vertex_set[f] = Set{Int}()
            for e in edges_in_face
                face_vertex_set[f] = union(face_vertex_set[f], edge_vertex_set[e])
            end
        end

        for c in 1:N_cells
            faces_in_cell = view(rowvals(E2T), nzrange(E2T, c))
            cell_vertex_set[c] = Set{Int}()
            for f in faces_in_cell
                cell_vertex_set[c] = union(cell_vertex_set[c], face_vertex_set[f])
            end
        end

        return new{DMT,SMT}(
            vertices, cell_centers,
            E0, E1, E2,
            E0T, E1T, E2T,
            edge_vertex_set,
            face_vertex_set,
            cell_vertex_set)
    end
end

n_vertices(complex::CellComplex) = size(complex.vertices, 2)
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
    edges_to_drop = ((complex.E0 * vertices_to_drop) .> 0) .|| vec(sum(complex.E0, dims=2) .< 2)
    faces_to_drop = (complex.E1 * edges_to_drop) .> 0
    volumes_to_drop = (complex.E2 * faces_to_drop) .> 0

    return SubComplex(complex, .!vertices_to_drop, .!edges_to_drop, .!faces_to_drop, .!volumes_to_drop)
end

function realize(sc::SubComplex{DMT,SMT,V}) where {DMT,SMT,V}
    E0 = sc.parent.E0[sc.edge_sub, sc.vertex_sub]
    E1 = sc.parent.E1[sc.face_sub, sc.edge_sub]
    E2 = sc.parent.E2[sc.volume_sub, sc.face_sub]
    vertices = sc.parent.vertices[:, sc.vertex_sub]
    cell_centers = sc.parent.cell_centers[:, sc.volume_sub]

    return CellComplex{DMT,SMT}(
        vertices, cell_centers,
        E0, E1, E2)
end

# o is opposite the desired normal direction
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