# Ω is a signed distance function
function bound_voronoi(unbounded_voronoi::CellComplex{DMT,SMT}, delaunay::CellComplex{DMT,SMT}, Ω::F; tol=1e-2) where {DMT,SMT,F}

    # Get subset of the voronoi with points outside Ω removed
    to_drop = .!WeightedCVT.in_Ω(unbounded_voronoi.vertices, Ω)
    pruned_voronoi = prune_points(unbounded_voronoi, to_drop)

    # Find edges that should exist based on faces present in the Delaunay tetrahedralization, but that do not exist in the pruned Voronoi diagram
    voronoi_edges_to_add = (((delaunay.E2T * pruned_voronoi.vertex_sub) .> 0) .&& .!pruned_voronoi.edge_sub)

    # There will also be just as many new points
    N_new_edges = sum(voronoi_edges_to_add)

    # Find starting points to form new edges
    # These are 
    i = 1
    new_points = DMT(zeros(3, N_new_edges))
    edge_indices = zeros(Int, N_new_edges)
    for e in 1:n_edges(unbounded_voronoi)
        if !voronoi_edges_to_add[e]
            continue
        end
        edge_indices[i] = e

        # Find the start vertex, the Voronoi vertex of the edges-to-add that is present in the subset
        vertices_in_edge = unbounded_voronoi.edge_vertex_set[e]
        to_drop = 0
        start_vertex::Int = 0
        for v in vertices_in_edge
            if pruned_voronoi.vertex_sub[v]
                start_vertex = v
            else
                to_drop = v
            end
        end

        # Get 
        # - the Delaunay face corresponding to the edge
        # - the Delaunay vertex opposite it
        #   - which is the tetrahedron vertex not in the face vertex set 
        # Use that to get the normal of the Delaunay face (direction vector for the edge)
        vertices_in_face = delaunay.face_vertex_set[e]
        vertices_in_tet = delaunay.cell_vertex_set[start_vertex]
        p0 = view(delaunay.vertices, :, collect(vertices_in_face)[1])
        p1 = view(delaunay.vertices, :, collect(vertices_in_face)[2])
        p2 = view(delaunay.vertices, :, collect(vertices_in_face)[3])
        po = view(delaunay.vertices, :, collect(setdiff(vertices_in_tet, vertices_in_face))[1])

        normal = triangle_normal(p0, p1, p2, po)
        init = SVector{3}(unbounded_voronoi.vertices[1, start_vertex],
            unbounded_voronoi.vertices[2, start_vertex],
            unbounded_voronoi.vertices[3, start_vertex])

        # Use normal vector and start vertex to get intersection with ∂Ω via SDF
        new_points[:, i] .= intersect_sdf(init, normal, Ω; tol=tol)

        # Need E0[e, i+Nv]=1, E0[e, start_vertex]=1
        # unbounded_voronoi.E0[e, i+n_vertices(unbounded_voronoi)] = 1 # this will be beyond the extant E0
        if to_drop > 0
            unbounded_voronoi.E0[e, to_drop] = 0
        end
        # Past this point, unbounded_voronoi should be considered potentially inconsistent

        i += 1
    end

    E0_new = hcat(dropzeros(unbounded_voronoi.E0), sparse(edge_indices, collect(1:N_new_edges), ones(N_new_edges), n_edges(unbounded_voronoi), N_new_edges))
    vertices_new::DMT = hcat(unbounded_voronoi.vertices, new_points)

    extant_indices = vcat(pruned_voronoi.vertex_sub, trues(N_new_edges))
    extant_edges = voronoi_edges_to_add .|| pruned_voronoi.edge_sub
    extant_faces = (unbounded_voronoi.E1 * extant_edges .> 0) .|| pruned_voronoi.face_sub


    # Close open face loops
    # Do this by finding which vertices are only incident to one face
    vertex_face_incidence = transpose(E0_new) * unbounded_voronoi.E1T
    dangling_vertices = ((Diagonal(extant_indices) * vertex_face_incidence) .== 1)

    # The faces which have such vertices incident to them will get a new edge connecting them
    open_faces = vec(sum(dangling_vertices, dims=1)) .> 0
    N_new_face_edges = sum(open_faces)
    new_face_cols = 1:N_new_face_edges
    extant_edges = vcat(extant_edges, trues(N_new_face_edges))
    E0_new = vcat(E0_new, sparse(repeat(1:N_new_face_edges, inner=2), rowvals(dangling_vertices), ones(2N_new_face_edges), N_new_face_edges, size(E0_new, 2)))
    E1_new = hcat(unbounded_voronoi.E1, sparse(findall(open_faces), new_face_cols, ones(N_new_face_edges), n_faces(unbounded_voronoi), N_new_face_edges))

    # Finally, draw some ray out from Delaunay boundary vertices (Voronoi generators) to ∂Ω, making a point p
    # Use p and those edges only incident to one face to close open cells with triangles
    dangling_edges = (Diagonal(extant_edges) * transpose(E1_new) * unbounded_voronoi.E2T) .== 1
    dangling_vertices = Diagonal(extant_indices) * transpose(E0_new) * dangling_edges
    open_cells = vec(sum(dangling_edges, dims=1)) .> 0
    N_new_surface_points = sum(open_cells)
    E0_newT = sparse(transpose(E0_new))

    i = 1
    new_surface_points = zeros(3, N_new_surface_points)
    new_edge_cols = Vector{Int}()
    N_new_surface_edges = 0

    new_face_cols = Vector{Int}()
    cell_row = Vector{Int}()

    # new_face_rows
    for c in 1:n_volumes(unbounded_voronoi)
        if !open_cells[c]
            continue
        end

        vertices_to_close = view(rowvals(dangling_vertices), nzrange(dangling_vertices, c))
        edges_to_close = view(rowvals(dangling_edges), nzrange(dangling_edges, c))
        # First, compute mean out towards surface direction 
        init = SVector{3}(unbounded_voronoi.cell_centers[1, c], unbounded_voronoi.cell_centers[2, c], unbounded_voronoi.cell_centers[3, c])
        dir = zero(SVector{3})
        for v::Int in vertices_to_close
            dir += SVector{3}(vertices_new[1, v], vertices_new[2, v], vertices_new[3, v])
        end
        dir = (dir / length(vertices_to_close)) - init
        dir /= norm(dir)

        # Find intersection point with ∂Ω
        new_surface_points[:, i] .= intersect_sdf(init, dir, Ω; tol=0.01)

        # Add edges between each vertex and this point, and close the triangular faces thus formed
        # for v in vertices_to_close
        inward_edge = Dict{Int,Int}()
        for e in edges_to_close
            vs = view(rowvals(E0_newT), nzrange(E0_newT, e))
            v1::Int = vs[1]
            v2::Int = vs[2]

            if !haskey(inward_edge, v1)
                N_new_surface_edges += 1
                push!(new_edge_cols, v1)
                push!(new_edge_cols, i + size(vertices_new, 2))
                inward_edge[v1] = N_new_surface_edges + size(E0_new, 1)
            end

            if !haskey(inward_edge, v2)
                N_new_surface_edges += 1
                push!(new_edge_cols, v2)
                push!(new_edge_cols, i + size(vertices_new, 2))
                inward_edge[v2] = N_new_surface_edges + size(E0_new, 1)
            end

            e1 = inward_edge[v1]
            e2 = inward_edge[v2]

            push!(new_face_cols, e)
            push!(new_face_cols, e1)
            push!(new_face_cols, e2)

            push!(cell_row, c)
        end

        i += 1
    end

    E0_closing = sparse(repeat(1:N_new_surface_edges, inner=2), new_edge_cols, ones(length(new_edge_cols)), N_new_surface_edges, N_new_surface_points + size(vertices_new, 2))
    E1_closing = sparse(repeat(1:N_new_surface_edges, inner=3), new_face_cols, ones(length(new_face_cols)), N_new_surface_edges, N_new_surface_edges + size(E1_new, 2))

    E0_new = vcat(
        hcat(E0_new, spzeros(size(E0_new, 1), N_new_surface_points)),
        E0_closing
    )

    E1_new = vcat(
        hcat(E1_new, spzeros(size(E1_new, 1), N_new_surface_edges)),
        E1_closing
    )

    E2_new = hcat(unbounded_voronoi.E2, sparse(cell_row, 1:N_new_surface_edges, ones(N_new_surface_edges), n_volumes(unbounded_voronoi), N_new_surface_edges))

    extant_indices = vcat(extant_indices, trues(N_new_surface_points))
    extant_edges = vcat(extant_edges, trues(N_new_surface_edges))
    extant_faces = vcat(extant_faces, trues(N_new_surface_edges))

    vertices_new = hcat(vertices_new, new_surface_points)

    return CellComplex{DMT,SMT}(
        vertices_new[:, extant_indices], unbounded_voronoi.cell_centers,
        E0_new[extant_edges, extant_indices], E1_new[extant_faces, extant_edges], E2_new[:, extant_faces]
    )
end

function bounded_voronoi(points::DMT, Ω::F; SMT=SparseMatrixCSC) where {DMT,F}
    extended_points = hcat([100 0 0 0; 0 100 0 0; 0 0 100 0], points)
    t = delaunay_tet(extended_points)

    r, dense_delaunay = condense_delaunay(t, extended_points; SMT=SMT)
    dense_voronoi = dual_complex(dense_delaunay)

    return bound_voronoi(dense_voronoi, dense_delaunay, Ω)
end