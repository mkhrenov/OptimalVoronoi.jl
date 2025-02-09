# Visualization

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
