# Incidence matrix representation for a 3-complex
struct CellComplex{DMT,SMT}
    vertices::DMT
    cell_centers::DMT

    E0::SMT # Vertex to edge incidence matrix
    E1::SMT # Edge to face incidence matrix
    E2::SMT # Face to volume incidence matrix
end

n_vertices(cell::CellComplex) = size(cell.vertices, 1)
n_edges(cell::CellComplex) = size(cell.E0, 1)
n_faces(cell::CellComplex) = size(cell.E1, 1)
n_volumes(cell::CellComplex) = size(cell.E2, 1)