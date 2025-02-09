function dual_complex(primal::CellComplex{DMT,SMT}) where {DMT,SMT}
    # Vertices of dual are circumcenters of cells/volumes of primal and vice-vers
    dual_vertices = copy(primal.cell_centers)
    dual_cell_centers = copy(primal.vertices)

    # Edges of dual are normal to faces of primal and connect incident volumes
    dual_E0 = SMT(transpose(primal.E2))

    # Faces of dual are normal to edges of primal and connect primal vertices (dual volumes)
    dual_E1 = SMT(transpose(primal.E1))

    # Volumes of dual are centered on vertices of the primal and have faces corresponding to each edge of the primal incident on primal vertex
    dual_E2 = SMT(transpose(primal.E0))

    return CellComplex{DMT,SMT}(
        dual_vertices, dual_cell_centers,
        dual_E0, dual_E1, dual_E2
    )
end