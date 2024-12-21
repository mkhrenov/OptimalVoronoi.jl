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

"""
    in_sphere(a, b, c, d, p)

Returns a positivie value if point `p` is within the sphere defined by points `a`, `b`, `c`, and `d`, negative otherwise.

For this to hold, `a`, `b`, `c`, and `d` must be ordered such that `orient(a, b, c, d)` returns a positive value.
"""
function in_sphere(a::VT, b::VT, c::VT, d::VT, p::VT) where {T,VT<:AbstractVector{T}}
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
