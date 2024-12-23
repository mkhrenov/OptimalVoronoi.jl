
struct DenseDelaunay{DIM,A1,A2}
    points::A1
    simplices::A2
    neighbors::A2

    function DenseDelaunay{DIM,A1,A2}(root_simplex::DelaunaySimplex{DIM}, points::A1) where {DIM,A1,A2}
        n_simplices = count_simplices(root_simplex)
        simplices = zeros(Int, DIM, n_simplices)
        neighbors = zeros(Int, DIM, n_simplices)

        simplex_map = collect_simplices!(simplices, root_simplex)
        collect_neighbors!(neighbors, root_simplex, simplex_map)

        simplices = A2(simplices)
        neighbors = A2(neighbors)

        new{DIM,A1,A2}(points, simplices, neighbors)
    end
end

function count_simplices(simplex::DelaunaySimplex{DIM}) where {DIM}
    scheduled = Set{DelaunaySimplex{DIM}}()
    tovisit = Stack{DelaunaySimplex{DIM}}()

    push!(tovisit, simplex)
    push!(scheduled, simplex)

    while !isempty(tovisit)
        current = pop!(tovisit)

        for neighbor in current.neighbors
            if neighbor ∉ scheduled
                push!(tovisit, neighbor)
                push!(scheduled, neighbor)
            end
        end
    end

    return length(scheduled)
end

function collect_simplices!(simplices, simplex::DelaunaySimplex{DIM}) where {DIM}
    scheduled = Set{DelaunaySimplex{DIM}}()
    tovisit = Stack{DelaunaySimplex{DIM}}()
    simplex_map = Dict{DelaunaySimplex{DIM},Int}()

    push!(tovisit, simplex)
    push!(scheduled, simplex)
    i = 1

    while !isempty(tovisit)
        current = pop!(tovisit)

        for neighbor in current.neighbors
            if neighbor ∉ scheduled
                push!(tovisit, neighbor)
                push!(scheduled, neighbor)
            end
        end

        for j in 1:DIM
            simplices[j, i] = current.vertices[j]
        end

        simplex_map[current] = i
        i += 1
    end

    return simplex_map
end

function collect_neighbors!(neighbors, simplex::DelaunaySimplex{DIM}, simplex_map) where {DIM}
    scheduled = Set{DelaunaySimplex{DIM}}()
    tovisit = Stack{DelaunaySimplex{DIM}}()

    push!(tovisit, simplex)
    push!(scheduled, simplex)

    while !isempty(tovisit)
        current = pop!(tovisit)

        for neighbor in current.neighbors
            if neighbor ∉ scheduled
                push!(tovisit, neighbor)
                push!(scheduled, neighbor)
            end
        end

        i = simplex_map[current]
        for j in 1:DIM
            neighbors[j, i] = simplex_map[current.neighbors[j]]
        end
    end

    return nothing
end