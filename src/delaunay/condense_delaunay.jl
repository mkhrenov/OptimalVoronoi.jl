function condense_delaunay(root_simplex::DelaunaySimplex{DIM}, vertices::DMT; SMT=SparseMatrixCSC{Int,Int}) where {DIM,DMT}
    simplex_list, simplex_map = enumerate_simplices(root_simplex)
    face_list, face_map, volumes_to_faces = enumerate_faces(simplex_list, simplex_map)
    edge_list, edge_map, faces_to_edges = enumerate_edges(face_list, face_map)

    reduced_vertices = vertices[:, 5:end]

    N_vertices = size(reduced_vertices, 2)
    N_edges = length(edge_list)
    N_faces = length(face_list)
    N_volumes = length(simplex_list)

    E2 = SMT(adjacency_list_to_incidence_matrix(volumes_to_faces, N_faces))
    E1 = SMT(adjacency_list_to_incidence_matrix(faces_to_edges, N_edges))
    E0 = SMT(adjacency_list_to_incidence_matrix(reshape(reinterpret(Int, edge_list), (2, :))', N_vertices))

    cell_centers = DMT(undef, 3, N_volumes)
    for (i, simplex) in enumerate(simplex_list)
        circumcenter!(view(cell_centers, :, i), simplex, reduced_vertices)
    end

    return CellComplex{DMT,SMT}(
        reduced_vertices, cell_centers,
        E0, E1, E2)
end

function enumerate_simplices(simplex::DelaunaySimplex{DIM}) where {DIM}
    scheduled = Set{DelaunaySimplex{DIM}}()
    tovisit = Stack{DelaunaySimplex{DIM}}()

    simplex_list = Vector{NTuple{4,Int}}()
    simplex_map = Dict{NTuple{4,Int},Int}()
    i = 0

    push!(tovisit, simplex)
    push!(scheduled, simplex)

    while !isempty(tovisit)
        current = pop!(tovisit)

        if !(1 in current.vertices || 2 in current.vertices || 3 in current.vertices || 4 in current.vertices)
            i += 1
            simp = tuple(
                current.vertices[1] - 4,
                current.vertices[2] - 4,
                current.vertices[3] - 4,
                current.vertices[4] - 4)
            simplex_map[simp] = i
            push!(simplex_list, simp)
        end

        for neighbor in current.neighbors
            if neighbor ∉ scheduled
                push!(tovisit, neighbor)
                push!(scheduled, neighbor)
            end
        end
    end

    return simplex_list, simplex_map
end


function enumerate_faces(simplex_list::Vector{NTuple{DIM,Int}}, simplex_map::Dict{NTuple{DIM,Int},Int}) where {DIM}
    face_list = Vector{NTuple{3,Int}}()
    face_map = Dict{NTuple{3,Int},Int}()
    volumes_to_faces = zeros(Int, length(simplex_list), 4)
    idx = 0

    temp_face_vec = zeros(Int, 3)

    for simplex in simplex_list
        simplex_idx = simplex_map[simplex]
        c = 1

        for i in 1:(DIM-2)
            for j in (i+1):(DIM-1)
                for k in (j+1):DIM
                    temp_face_vec[1] = simplex[i]
                    temp_face_vec[2] = simplex[j]
                    temp_face_vec[3] = simplex[k]

                    sort!(temp_face_vec)
                    face = tuple(temp_face_vec[1], temp_face_vec[2], temp_face_vec[3])::NTuple{3,Int}

                    if !haskey(face_map, face)
                        idx += 1
                        face_map[face] = idx
                        push!(face_list, face)
                    end

                    face_idx = face_map[face]
                    volumes_to_faces[simplex_idx, c] = face_idx
                    c += 1
                end
            end
        end
    end

    return face_list, face_map, volumes_to_faces
end

function enumerate_edges(face_list::Vector{NTuple{DIM,Int}}, face_map::Dict{NTuple{DIM,Int},Int}) where {DIM}
    edge_list = Vector{NTuple{2,Int}}()
    edge_map = Dict{NTuple{2,Int},Int}()
    faces_to_edges = zeros(Int, length(face_list), 3)
    idx = 0

    temp_edge_vec = zeros(Int, 2)

    for face in face_list
        face_idx = face_map[face]
        c = 1

        for i in 1:(DIM-1)
            for j in (i+1):DIM
                temp_edge_vec[1] = face[i]
                temp_edge_vec[2] = face[j]

                sort!(temp_edge_vec)
                edge = NTuple{2,Int}(temp_edge_vec)

                if !haskey(edge_map, edge)
                    idx += 1
                    edge_map[edge] = idx
                    push!(edge_list, edge)
                end

                edge_idx = edge_map[edge]
                faces_to_edges[face_idx, c] = edge_idx
                c += 1
            end
        end
    end

    return edge_list, edge_map, faces_to_edges
end

function adjacency_list_to_incidence_matrix(adj_list, n_cols)
    n_rows, k = size(adj_list)

    rv = zeros(Int, n_rows * k)
    cv = zeros(Int, n_rows * k)
    vv = ones(Int, n_rows * k)

    for i in 1:n_rows
        for j in 1:k
            r = i
            c = adj_list[i, j]

            rv[(i-1)*k+j] = r
            cv[(i-1)*k+j] = c
        end
    end

    return sparse(rv, cv, vv, n_rows, n_cols)
end


function circumcenter!(center_vec, simplex::NTuple{4,Int}, vertices)
    v1 = SVector{3}(view(vertices, :, simplex[1]))
    v2 = SVector{3}(view(vertices, :, simplex[2]))
    v3 = SVector{3}(view(vertices, :, simplex[3]))
    v4 = SVector{3}(view(vertices, :, simplex[4]))
    a = det(SMatrix{4,4}(
        v1[1], v1[2], v1[3], 1,
        v2[1], v2[2], v2[3], 1,
        v3[1], v3[2], v3[3], 1,
        v4[1], v4[2], v4[3], 1,
    ))

    r1 = v1[1]^2 + v1[2]^2 + v1[3]^2
    r2 = v2[1]^2 + v2[2]^2 + v2[3]^2
    r3 = v3[1]^2 + v3[2]^2 + v3[3]^2
    r4 = v4[1]^2 + v4[2]^2 + v4[3]^2

    Dx = det(SMatrix{4,4}(
        r1, v1[2], v1[3], 1,
        r2, v2[2], v2[3], 1,
        r3, v3[2], v3[3], 1,
        r4, v4[2], v4[3], 1,
    ))

    Dy = -det(SMatrix{4,4}(
        r1, v1[1], v1[3], 1,
        r2, v2[1], v2[3], 1,
        r3, v3[1], v3[3], 1,
        r4, v4[1], v4[3], 1,
    ))

    Dz = det(SMatrix{4,4}(
        r1, v1[1], v1[2], 1,
        r2, v2[1], v2[2], 1,
        r3, v3[1], v3[2], 1,
        r4, v4[1], v4[2], 1,
    ))

    c = det(SMatrix{4,4}(
        r1, v1[1], v1[2], v1[3],
        r2, v2[1], v2[2], v2[3],
        r3, v3[1], v3[2], v3[3],
        r4, v4[1], v4[2], v4[3],
    ))

    center_vec[1] = Dx / (2 * a)
    center_vec[2] = Dy / (2 * a)
    center_vec[3] = Dz / (2 * a)

    return √(Dx^2 + Dy^2 + Dz^2 - 4 * a * c) / abs(2 * a)
end