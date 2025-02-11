# Ω is a signed distance function
function bound_voronoi(unbounded_voronoi::CellComplex{DMT,SMT}, delaunay::CellComplex{DMT,SMT}, Ω::F) where {DMT,SMT,F}

    # Get subset of the voronoi with points outside Ω removed
    to_drop = .!WeightedCVT.in_Ω(unbounded_voronoi.vertices, Ω)
    pruned_voronoi = prune_points(unbounded_voronoi, to_drop)

    # Find edges that should exist based on faces present in the Delaunay tetrahedralization, but that do not exist in the pruned Voronoi diagram
    voronoi_edges_to_add = (((delaunay.E2T * pruned_voronoi.vertex_sub) .> 0) .&& .!pruned_voronoi.edge_sub)

    # There will also be just as many new points
    N_new_edges = sum(voronoi_edges_to_add)

    # Find normals and starting points to form new edges
    normals = zeros(3, N_new_edges)
    for i in 1:N_new_edges
        
    end
    start_points = 0

    # Find termination points of new edges through search on SDF, drawing rays normal to 
    # Delaunay boundary faces from tetrahedron circumcenters (Voronoi vertices) until they intersect ∂Ω at a point
    new_points = 0

    # Create entries of E0 to add new edges

    # Figure out what edges should be added on ∂Ω to connect new points, and close loops to form faces between Voronoi cells
    # (Brute force geometric -- the coplanar ones)

    # Finally, draw some ray out from Delaunay boundary vertices (Voronoi generators) to ∂Ω,
    # and use it and preceding edges to form triangular faces
end