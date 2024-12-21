using LinearAlgebra
using StaticArrays
using DataStructures
using GLMakie

##### Data Structures #####

struct Delaunay{PAT,SAT,DIM}
    points::PAT
    simplices::SAT
    neighbors::SAT
end

struct DelaunaySimplex{DIM}
    vertices::NTuple{DIM,Int}
    neighbors::Vector{DelaunaySimplex{DIM}}

    function DelaunaySimplex{4}(v1::Int, v2::Int, v3::Int, v4::Int)
        neighbors = Vector{DelaunaySimplex{4}}(undef, 4)

        simplex = new{4}((v1, v2, v3, v4), neighbors)
        simplex.neighbors[1] = simplex
        simplex.neighbors[2] = simplex
        simplex.neighbors[3] = simplex
        simplex.neighbors[4] = simplex
    end
end

const DelaunayTet = DelaunaySimplex{4}

include("predicates.jl")
include("flips.jl")

function replace_neighbor!(simplex::DelaunaySimplex{DIM}, old_neighbor::DelaunaySimplex{DIM}, new_neighbor::DelaunaySimplex{DIM}) where {DIM}
    for i in 1:DIM
        if simplex.neighbors[i] === old_neighbor
            simplex.neighbors[i] = new_neighbor
            break
        end
    end
end

##### Point Location #####

function walk(points::PMT, σ::DelaunaySimplex{DIM}, p::Int) where {DIM,T<:Number,PMT<:AbstractMatrix{T}}

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
    v1 = 0
    v2 = 0
    v3 = 0
    q = 0

    # Find the points that form the interface between n and σ, and the point `q` in σ which is opposite the interface
    if n == 1
        v1 = σ.vertices[1]
        v2 = σ.vertices[3]
        v3 = σ.vertices[2]
        q = σ.vertices[4]
    elseif n == 2
        v1 = σ.vertices[1]
        v2 = σ.vertices[3]
        v3 = σ.vertices[2]
        q = σ.vertices[4]
    elseif n == 3
        v1 = σ.vertices[1]
        v2 = σ.vertices[3]
        v3 = σ.vertices[2]
        q = σ.vertices[4]
    elseif n == 4
        v1 = σ.vertices[1]
        v2 = σ.vertices[3]
        v3 = σ.vertices[2]
        q = σ.vertices[4]
    end

    # Get the orientation of `q` wrt the interface
    oq = orient(@view(points[:, v1]), @view(points[:, v2]), @view(points[:, v3]), @view(points[:, q]))

    # Get the orientation of `p` wrt the interface
    op = orient(@view(points[:, v1]), @view(points[:, v2]), @view(points[:, v3]), @view(points[:, p]))

    # If they have opposite signs, n is closer to p than σ
    return oq * op < 0
end

##### Delaunay Triangulation #####

function insert_one_point(DT::Delaunay, p::Int)

end

function delaunay()

end



function viz(simplex::DelaunaySimplex{DIM}, points) where {DIM}
    visited = Set{DelaunaySimplex{DIM}}()
    tovisit = Stack{DelaunaySimplex{DIM}}()

    fig, ax, plot = scatter(points)

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
end

function viz!(simplex::DelaunaySimplex{DIM}, points) where {DIM}
    segment = zeros(DIM-1, 2)

    for i in 1:(DIM-1)
        segment[:, 1] .= @view points[:, simplex.vertices[i]]

        for j in (i+1):DIM
            segment[:, 2] .= @view points[:, simplex.vertices[j]]
            lines!(segment, color=:red)
        end
    end
end


function viz(DT::Delaunay{PAT,SAT,DIM}) where {PAT,SAT,DIM}
    fig, ax, plt = scatter(DT.points)
    n_simplices = size(DT.simplices, 2)
    segment = zeros(3, 2)

    for k in 1:n_simplices
        for i in 1:DIM
            segment[:, 1] .= @view DT.points[:, DT.simplices[k, i]]

            for j in (i+1):(DIM+1)
                segment[:, 2] .= @view DT.points[:, DT.simplices[k, j]]
                lines!(segment, color=:red)
            end
        end
    end

    display(fig)
end