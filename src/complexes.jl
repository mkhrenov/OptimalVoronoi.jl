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

n_verts(complex::CellComplex) = size(complex.vertices, 2)
n_edges(complex::CellComplex) = size(complex.E0, 1)
n_faces(complex::CellComplex) = size(complex.E1, 1)
n_cells(complex::CellComplex) = size(complex.E2, 1)

# Iterators for CSC matrices

function verts_of_edge(complex::CellComplex{DMT,SparseMatrixCSC{Int,Int}}, edge::Int) where {DMT}
    return view(rowvals(complex.E0T), nzrange(complex.E0T, edge))
end

function edges_of_face(complex::CellComplex{DMT,SparseMatrixCSC{Int,Int}}, face::Int) where {DMT}
    return view(rowvals(complex.E1T), nzrange(complex.E1T, face))
end

function faces_of_cell(complex::CellComplex{DMT,SparseMatrixCSC{Int,Int}}, cell::Int) where {DMT}
    return view(rowvals(complex.E2T), nzrange(complex.E2T, cell))
end

function edges_of_vert(complex::CellComplex{DMT,SparseMatrixCSC{Int,Int}}, vert::Int) where {DMT}
    return view(rowvals(complex.E0), nzrange(complex.E0, vert))
end

function faces_of_edge(complex::CellComplex{DMT,SparseMatrixCSC{Int,Int}}, edge::Int) where {DMT}
    return view(rowvals(complex.E1), nzrange(complex.E1, edge))
end

function cells_of_face(complex::CellComplex{DMT,SparseMatrixCSC{Int,Int}}, face::Int) where {DMT}
    return view(rowvals(complex.E2), nzrange(complex.E2, face))
end

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
    volume_sub = falses(n_cells(complex))
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

    n = v1 × v2

    if (ov - p0) ⋅ n > 0.0
        n = n * -1
    end

    return n / norm(n)
end