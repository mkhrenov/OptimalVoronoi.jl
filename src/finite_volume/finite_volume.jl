function mesh_fv_matrix_vector(complex::CellComplex{DMT,SMT}, ρ, k, cₚ, h) where {DMT,SMT}
    α = k / ρ / cₚ
    volumes = complex_volumes(complex)

    Adj, e = mesh_adjacency_matrix_vector(complex)

    Deg = Diagonal(vec(sum(Adj, dims=2)))
    L = Deg - Adj
    # Conduction given by -L ⋅ k / (ρ ⋅ cₚ ⋅ Vᵢ)
    Aᵢₙ = -α * L ./ volumes

    e ./= volumes
    e .*= (h / ρ / cₚ)
    Aₑₓₜ = Diagonal(e)
    A = Aᵢₙ - Aₑₓₜ

    return A, e
end

function mesh_adjacency_matrix_vector(complex::CellComplex{DMT,SMT}) where {DMT,SMT}
    # Aᵢⱼ / lᵢⱼ, j ∈ Nᵢ
    rows::Vector{Int} = []
    cols::Vector{Int} = []
    vals::Vector{Float64} = []

    e = zeros(n_cells(complex))

    for i in 1:n_cells(complex)
        for face in faces_of_cell(complex, i)
            j = other_cell_of_face(complex, face, i)

            if j == 0
                e[i] += face_area(complex, face)
            else
                push!(rows, i)
                push!(cols, j)
                push!(vals, face_area(complex, face) / cell_separation(complex, i, j))
            end
        end
    end

    Adj = sparse(rows, cols, vals)
    return Adj, e
end

function mesh_fv_matrix_vector(points::MT, domain, Ω::G, l, ρ, k, cₚ, h) where {MT,G}
    N_cells = size(points, 2)

    α = k / ρ / cₚ

    volumes = MT(zeros(1, N_cells))
    cell_volume_integrals!(volumes, (p, i) -> 1, domain)
    volumes = vec(volumes)
    volumes .*= l^3

    Adj, e = mesh_adjacency_matrix_vector(points, domain, Ω, l)
    Deg = Diagonal(vec(sum(Adj, dims=2)))
    L = typeof(Adj)(Deg) .- Adj

    # Conduction given by -L ⋅ k / (ρ ⋅ cₚ ⋅ Vᵢ)
    Aᵢₙ = -α * Diagonal(1 ./ volumes) * L

    e ./= volumes
    e .*= (h / ρ / cₚ)
    Aₑₓₜ = Diagonal(e)
    A = Aᵢₙ - Aₑₓₜ

    return A, e
end

function mesh_adjacency_matrix_vector(points::MT, domain, Ω::G, l) where {MT,G}
    N_cells = size(points, 2)
    A = MT(zeros(N_cells, N_cells))
    e = vec(MT(zeros(1, N_cells)))

    color_voronoi!(domain, points)
    adjacency_matrix_vector!(A, e, domain)

    Adj = sparse(A)
    neighbor_surface_integrals!(Adj, (p, i, j) -> l / point_separation(points, i, j), domain, points)

    exterior_surface_integrals!(e, (p) -> l^2, Ω, domain, points)

    return (Adj + transpose(Adj)) / 2, e
end

function point_separation(points, i, j)
    p_i = SVector{3}(points[1, i], points[2, i], points[2, i])
    p_j = SVector{3}(points[1, j], points[2, j], points[2, j])
    return norm(p_i - p_j)
end