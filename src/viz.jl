# Visualization

function viz(complex::CellComplex{DMT,SMT}; edge_color=:green, cell_colors=nothing) where {DMT,SMT}
    fig = Figure()
    lscene = LScene(fig[1, 1])
    viz!(complex; edge_color=edge_color, cell_colors=cell_colors)
end

function viz!(complex::CellComplex{DMT,SMT}; edge_color=:green, cell_colors=nothing) where {DMT,SMT}
    scatter!(complex.vertices)
    scatter!(complex.cell_centers)

    if isnothing(cell_colors)
        cell_colors = [RGBf(rand(3)...) for _ in 1:n_cells(complex)]
    else
        cell_colors = floor.(Int, 256 .* (cell_colors .- minimum(cell_colors)) ./ (maximum(cell_colors) + 0.01 - minimum(cell_colors))) .+ 1
        cell_colors = cgrad(:default; alpha=0.5)[cell_colors]
    end

    segment = zeros(3, 2)
    rowvec = rowvals(complex.E0T)

    for i in 1:n_edges(complex)
        if length(nzrange(complex.E0T, i)) != 2
            continue
        end

        for (c, j) in enumerate(nzrange(complex.E0T, i))
            v = rowvec[j]
            segment[:, c] .= @view complex.vertices[:, v]
        end

        lines!(segment, color=edge_color)
    end

    triangle = zeros(3, 3)
    triangle_faces = [1 2 3]
    for face in 1:n_faces(complex)
        edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, face))
        v0 = view(rowvals(complex.E0T), nzrange(complex.E0T, edges_in_face[1]))[1]
        triangle[:, 1] .= @view complex.vertices[:, v0]

        face_color = zero(RGBf)

        cells_of_face = view(rowvals(complex.E2), nzrange(complex.E2, face))
        for cell in cells_of_face
            face_color += cell_colors[cell] / length(cells_of_face)
        end

        for edge in @view edges_in_face[2:end]
            vertices_in_edge = view(rowvals(complex.E0T), nzrange(complex.E0T, edge))

            if length(vertices_in_edge) != 2
                continue
            end

            for (c, vertex) in enumerate(vertices_in_edge)
                triangle[:, c+1] .= @view complex.vertices[:, vertex]
            end

            mesh!(triangle, triangle_faces, color=face_color, alpha=0.5, backlight=1.0)
        end
    end

    display(current_figure())
end

function viz(subcomplex::SubComplex{DMT,SMT,V}; edge_color=:green) where {DMT,SMT,V}
    fig = Figure()
    lscene = LScene(fig[1, 1])
    viz!(subcomplex; edge_color=edge_color)
end

function viz!(subcomplex::SubComplex{DMT,SMT,V}; edge_color=:green) where {DMT,SMT,V}
    complex = subcomplex.parent

    scatter!(complex.vertices[:, subcomplex.vertex_sub])
    scatter!(complex.cell_centers[:, subcomplex.volume_sub])

    segment = zeros(3, 2)
    rowvec = rowvals(complex.E0T)

    for i in 1:n_edges(complex)
        if !subcomplex.edge_sub[i] || length(nzrange(complex.E0T, i)) != 2
            continue
        end

        for (c, j) in enumerate(nzrange(complex.E0T, i))
            v = rowvec[j]
            segment[:, c] .= @view complex.vertices[:, v]
        end

        lines!(segment, color=edge_color)
    end

    triangle = zeros(3, 3)
    triangle_faces = [1 2 3]
    for i in 1:n_faces(complex)
        if !subcomplex.face_sub[i]
            continue
        end

        edges_in_face = view(rowvals(complex.E1T), nzrange(complex.E1T, i))
        v0 = view(rowvals(complex.E0T), nzrange(complex.E0T, edges_in_face[1]))[1]
        triangle[:, 1] .= @view complex.vertices[:, v0]

        face_color = RGBf(rand(3)...)

        for j in @view edges_in_face[2:end]
            vertices_in_edge = view(rowvals(complex.E0T), nzrange(complex.E0T, j))

            if length(vertices_in_edge) != 2
                continue
            end

            for (c, j) in enumerate(vertices_in_edge)
                triangle[:, c+1] .= @view complex.vertices[:, j]
            end

            mesh!(triangle, triangle_faces, color=face_color, alpha=0.5, backlight=1.0)
        end
    end

    display(current_figure())
end


function viz(simplex::DelaunaySimplex{DIM}, points) where {DIM}
    visited = Set{DelaunaySimplex{DIM}}()
    tovisit = Stack{DelaunaySimplex{DIM}}()

    fig, ax, plot = scatter(@view(points[:, 5:end]))

    push!(tovisit, simplex)

    while !isempty(tovisit)
        current = pop!(tovisit)
        push!(visited, current)

        viz!(current, points)

        for neighbor in current.neighbors
            if neighbor ∉ visited
                push!(tovisit, neighbor)
            end
        end
    end

    display(fig)
    return fig, ax, plot
end

function viz!(simplex::DelaunaySimplex{DIM}, points; color=:red) where {DIM}
    segment = zeros(DIM - 1, 2)

    for i in 1:(DIM-1)
        if simplex.vertices[i] ≤ 4
            continue
        end

        segment[:, 1] .= @view points[:, simplex.vertices[i]]

        for j in (i+1):DIM
            if simplex.vertices[j] ≤ 4
                continue
            end

            segment[:, 2] .= @view points[:, simplex.vertices[j]]
            lines!(segment, color=color)
        end
    end
end
