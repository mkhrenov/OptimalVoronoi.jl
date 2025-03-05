##### Data Structures #####

struct DelaunaySimplex{DIM}
    vertices::NTuple{DIM,Int}
    neighbors::Vector{DelaunaySimplex{DIM}}

    function DelaunaySimplex{4}(v1::Int, v2::Int, v3::Int, v4::Int)
        neighbors = Vector{DelaunaySimplex{4}}(undef, 4)

        simplex::DelaunaySimplex{4} = new{4}(tuple(v1, v2, v3, v4), neighbors)
        simplex.neighbors[1] = simplex
        simplex.neighbors[2] = simplex
        simplex.neighbors[3] = simplex
        simplex.neighbors[4] = simplex
    end
end

const DelaunayTet = DelaunaySimplex{4}

Base.show(io::IO, s::DelaunaySimplex) = print(io, "$(s.vertices)")

function simplex_is(s::DelaunaySimplex{DIM}, ps::Vararg{Int,DIM}) where {DIM}
    for p in ps
        if p ∉ s.vertices
            return false
        end
    end

    return true
end

function replace_neighbor!(simplex::DelaunaySimplex{DIM}, old_neighbor::DelaunaySimplex{DIM}, new_neighbor::DelaunaySimplex{DIM}) where {DIM}
    for i in 1:DIM
        if simplex.neighbors[i] === old_neighbor
            simplex.neighbors[i] = new_neighbor
            break
        end
    end
end

include("predicates.jl")
include("flips.jl")
include("condense_delaunay.jl")

##### Point Location #####

function walk(σ::DelaunaySimplex{DIM}, p::Int, points) where {DIM}
    converged = false
    while !converged
        converged = true

        for (n, neighbor) in enumerate(σ.neighbors)
            if neighbor == σ
                continue
            end

            if closertothan(p, n, σ, points)
                σ = neighbor
                converged = false
                break
            end
        end
    end

    return σ
end

function closertothan(p::Int, n::Int, σ::DelaunayTet, points::PMT) where {T<:Number,PMT<:AbstractMatrix{T}}
    # Find the points that form the interface between n and σ, and the point `q` in σ which is opposite the interface
    v1, v2, v3, q = get_commonface_and_opposite(σ, n)

    # Get the orientation of `q` wrt the interface
    oq = orient(@view(points[:, v1]), @view(points[:, v2]), @view(points[:, v3]), @view(points[:, q]))

    # Get the orientation of `p` wrt the interface
    op = orient(@view(points[:, v1]), @view(points[:, v2]), @view(points[:, v3]), @view(points[:, p]))

    # If they have opposite signs, n is closer to p than σ
    return oq * op < 0
end


##### Delaunay Tetrahedralization #####

function insert_one_point(genesis::DelaunayTet, p::Int, points, stack::Stack{DelaunayTet})
    τ = walk(genesis, p, points)
    t1, t2, t3, t4 = flip14(τ, p, points)

    deleted = Set{DelaunayTet}()
    push!(stack, t1)
    push!(stack, t2)
    push!(stack, t3)
    push!(stack, t4)

    while !isempty(stack)
        τ = pop!(stack)
        τₐ = get_neighbor_opposite(τ, p)

        if τ === τₐ || τ ∈ deleted || τₐ ∈ deleted
            continue
        end

        a, b, c, d = get_commonface_and_opposite(τₐ, τ)

        if in_circumsphere(τ, d, points)
            # Case 1?
            if line_intersects_triangle(p, d, a, b, c, points)
                t1, t2, t3 = flip23(τ, τₐ, points)
                push!(deleted, τ)
                push!(deleted, τₐ)

                push!(stack, t1)
                push!(stack, t2)
                push!(stack, t3)
            else#if !line_intersects_triangle(p, d, a, b, c, points)
                # check if tet opposite c is pdab
                pdab = get_neighbor_opposite(τ, c)
                if simplex_is(pdab, p, d, a, b)
                    t1, t2 = flip32(τ, τₐ, pdab, points)
                    push!(deleted, τ)
                    push!(deleted, τₐ)
                    push!(deleted, pdab)

                    push!(stack, t1)
                    push!(stack, t2)
                    continue
                end

                # check if tet opposite b is pdca
                pdca = get_neighbor_opposite(τ, b)
                if simplex_is(pdca, p, d, c, a)
                    t1, t2 = flip32(τ, τₐ, pdca, points)
                    push!(deleted, τ)
                    push!(deleted, τₐ)
                    push!(deleted, pdca)

                    push!(stack, t1)
                    push!(stack, t2)
                    continue
                end

                # check if tet opposite a is pdbc
                pdbc = get_neighbor_opposite(τ, a)
                if simplex_is(pdbc, p, d, b, c)
                    t1, t2 = flip32(τ, τₐ, pdbc, points)
                    push!(deleted, τ)
                    push!(deleted, τₐ)
                    push!(deleted, pdbc)

                    push!(stack, t1)
                    push!(stack, t2)
                    continue
                end
            end
        end
    end

    return t1
end

function delaunay_tet(points)
    DIM = size(points, 1)
    N = size(points, 2)
    @assert DIM == 3

    t = DelaunaySimplex{4}(1, 2, 3, 4)
    stack = Stack{DelaunayTet}()

    for i in (DIM+2):N
        t = insert_one_point(t, i, points, stack)
    end

    return t
end

function get_neighbor_opposite(t::DelaunayTet, p::Int)
    for i in 1:4
        if t.vertices[i] == p
            return t.neighbors[5-i]
        end
    end

    return t
end

function get_commonface_and_opposite(t::DelaunayTet, n::DelaunayTet)
    nn = get_neighbor_number(t, n)
    return get_commonface_and_opposite(t, nn)
end

function get_commonface_and_opposite(t::DelaunayTet, n::Int)
    if n == 1
        v1 = t.vertices[1]
        v2 = t.vertices[3]
        v3 = t.vertices[2]
        q = t.vertices[4]
        return v1, v2, v3, q
    elseif n == 2
        v1 = t.vertices[1]
        v2 = t.vertices[2]
        v3 = t.vertices[4]
        q = t.vertices[3]
        return v1, v2, v3, q
    elseif n == 3
        v1 = t.vertices[1]
        v2 = t.vertices[4]
        v3 = t.vertices[3]
        q = t.vertices[2]
        return v1, v2, v3, q
    else
        v1 = t.vertices[2]
        v2 = t.vertices[3]
        v3 = t.vertices[4]
        q = t.vertices[1]
        return v1, v2, v3, q
    end
end

function get_neighbor_number(t::DelaunayTet, n::DelaunayTet)
    for i in 1:4
        if t.neighbors[i] === n
            return i
        end
    end

    return 0
end

function is_delaunay(simplex::DelaunaySimplex{DIM}, points) where {DIM}
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