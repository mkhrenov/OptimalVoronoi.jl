using LinearAlgebra
using StaticArrays

include("predicates.jl")
include("flips.jl")


##### Data Structures #####

struct Delaunay{PAT,SAT,DIM}
    n_simplices::Int
    n_points::Int

    points::PAT
    simplices::SAT
    neighbors::SAT
end


##### Point Location #####

function walk(DT::Delaunay{PAT,SAT,DIM}, p::VT) where {T,VT<:AbstractVector{T},PAT,SAT,DIM}
    neighbors = DT.neighbors
    σ = rand(1:DT.n_simplices)

    converged = false
    while !converged
        converged = true

        for neighbor in @view(neighbors[σ, :])
            if closertothan(p, neighbor, σ, DT)
                σ = neighbor
                converged = false
                break
            end
        end
    end

    return σ
end

function closertothan(p::Int, n::Int, σ::Int, DT::Delaunay{PAT,SAT,DIM}) where {T,VT<:AbstractVector{T},PAT,SAT,DIM}
    sps = DT.simplices
    pts = DT.points

    c = 0
    v1 = 0
    v2 = 0
    v3 = 0
    q = 0

    # Find the points that form the interface between n and σ, and the point `q` in σ which is opposite the interface
    for i in @view(sps[σ, :])
        odd_one_out = true

        for j in @view(sps[n, :])
            if i == j
                if c == 0
                    v1 = i
                    c += 1
                    odd_one_out = false
                    break
                elseif c == 1
                    v2 = i
                    c += 1
                    odd_one_out = false
                    break
                elseif c == 2
                    v3 = i
                    c += 1
                    odd_one_out = false
                    break
                end
            end
        end

        if odd_one_out
            q = i
        end
    end

    # Get the orientation of `q` wrt the interface
    oq = orient(@view(pts[v1, :]), @view(pts[v2, :]), @view(pts[v3, :]), @view(pts[q, :]))

    # Get the orientation of `p` wrt the interface
    op = orient(@view(pts[v1, :]), @view(pts[v2, :]), @view(pts[v3, :]), @view(pts[p, :]))

    # If they have opposite signs, n is closer to p than σ
    return oq * op < 0
end

##### Delaunay Triangulation #####

function insert_one_point(DT::Delaunay, p::VT) where {T,VT<:AbstractVector{T}}

end

function delaunay()

end
