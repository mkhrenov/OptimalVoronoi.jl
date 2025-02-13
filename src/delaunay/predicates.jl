##### Predicates #####

"""
    orient(a, b, c, p)

Returns a positive value if point `p` is above the plane defined by points `a`, `b`, and `c`, negative otherwise.
"""
function orient(a::VT, b::VT, c::VT, p::VT) where {T,VT<:AbstractVector{T}}
    return det(SMatrix{4,4,T}(
        a[1], a[2], a[3], T(1.0),
        b[1], b[2], b[3], T(1.0),
        c[1], c[2], c[3], T(1.0),
        p[1], p[2], p[3], T(1.0),
    ))
end

function orient(a::Int, b::Int, c::Int, p::Int, points)
    return orient(
        @view(points[:, a]),
        @view(points[:, b]),
        @view(points[:, c]),
        @view(points[:, p])
    )
end

"""
    in_sphere(a, b, c, d, p)

Returns a positivie value if point `p` is within the sphere defined by points `a`, `b`, `c`, and `d`, negative otherwise.

For this to hold, `a`, `b`, `c`, and `d` must be ordered such that `orient(a, b, c, d)` returns a positive value.
"""
function in_sphere(a::VT, b::VT, c::VT, d::VT, p::VT) where {T,VT<:AbstractVector{T}}
    @assert orient(a, b, c, d) > 0.0
    return det(
        SMatrix{5,5,T}(
            a[1], a[2], a[3], a[1]^2 + a[2]^2 + a[3]^2, T(1.0),
            b[1], b[2], b[3], b[1]^2 + b[2]^2 + b[3]^2, T(1.0),
            c[1], c[2], c[3], c[1]^2 + c[2]^2 + c[3]^2, T(1.0),
            d[1], d[2], d[3], d[1]^2 + d[2]^2 + d[3]^2, T(1.0),
            p[1], p[2], p[3], p[1]^2 + p[2]^2 + p[3]^2, T(1.0),
        ),
    )
end


function in_circumsphere(t::DelaunayTet, p::Int, points)
    return in_sphere(
        @view(points[:, t.vertices[1]]),
        @view(points[:, t.vertices[2]]),
        @view(points[:, t.vertices[3]]),
        @view(points[:, t.vertices[4]]),
        @view(points[:, p])
    ) > 0.0
end


function line_intersects_triangle(p1::VT, p2::VT, a::VT, b::VT, c::VT) where {T,VT<:AbstractVector{T}}
    p1v = SVector{3,T}(p1[1], p1[2], p1[3])
    p2v = SVector{3,T}(p2[1], p2[2], p2[3])
    av = SVector{3,T}(a[1], a[2], a[3])
    bv = SVector{3,T}(b[1], b[2], b[3])
    cv = SVector{3,T}(c[1], c[2], c[3])

    dir = p2v - p1v
    origin = p1v

    e1 = bv - av
    e2 = cv - av

    n = e1 × e2
    det = -dot(dir, n)
    invdet = 1.0 / det

    a0 = origin - av
    da0 = a0 × dir

    u = dot(e2, da0) * invdet
    v = -dot(e1, da0) * invdet
    t = dot(a0, n) * invdet

    return (abs(det) ≥ 1e-6) && (t ≥ 0.0) && (u ≥ 0.0) && (v ≥ 0.0) && (u + v ≤ 1.0)
end

function line_intersects_triangle(p1::Int, p2::Int, a::Int, b::Int, c::Int, points)
    return line_intersects_triangle(
        @view(points[:, p1]),
        @view(points[:, p2]),
        @view(points[:, a]),
        @view(points[:, b]),
        @view(points[:, c])
    )
end
