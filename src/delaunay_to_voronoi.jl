
function count_cells(simplex::DelaunaySimplex{DIM}, points) where {DIM}
    scheduled = Set{DelaunaySimplex{DIM}}()
    tovisit = Stack{DelaunaySimplex{DIM}}()

    push!(tovisit, simplex)
    push!(scheduled, simplex)

    while !isempty(tovisit)
        current = pop!(tovisit)

        for i in 1:size(points, 2)
            if i ∉ current.vertices && in_circumsphere(current, i, points)
                return false
            end
        end

        for neighbor in current.neighbors
            if neighbor ∉ scheduled
                push!(tovisit, neighbor)
                push!(scheduled, neighbor)
            end
        end
    end

    return true
end