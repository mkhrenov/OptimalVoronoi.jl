function mesh_fv_matrix_vector(complex::CellComplex{DMT,SMT}, ρ, k, cₚ, h) where {DMT,SMT}
    α = k / ρ / cₚ
    volumes = complex_volumes(complex)

    Adj, e = mesh_adjacency_matrix_vector(complex)

    Deg = Diagonal(vec(sum(Adj, dims=2)))
    L = Deg - Adj
    # Conduction given by -L ⋅ k / (ρ ⋅ cₚ ⋅ Vᵢ)
    Aᵢₙ = -α * L ./ volumes

    e ./= volumes
    e .*= (h / ρ)
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