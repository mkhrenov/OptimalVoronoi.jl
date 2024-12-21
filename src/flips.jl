
# ##### Flips #####

# """
#     flip(tet1, tet2)


# """
# function flip(tet1::Tetrahedron, tet2::Tetrahedron)
#     if case1
#         flip23(tet1, tet2)
#     elseif case2
#         flip32(tet1, tet2)
#     elseif case3
#         flip44()
#     elseif case4
#         flip23()
#     end
# end

# function flip23(tet1::Tetrahedron, tet2::Tetrahedron)
#     a, b, c, d = tet1.vertices
#     b, c, d, e = tet2.vertices
#     # TODO: Figure out how to match orientations

#     return Tetrahedron(a, c, d, e), Tetrahedron(a, d, b, e), Tetrahedron(a, b, c, e)
# end

# function flip32(tet1::Tetrahedron, tet2::Tetrahedron, tet3::Tetrahedron)
#     a, c, d, e = tet1.vertices
#     a, d, b, e = tet2.vertices
#     a, b, c, e = tet3.vertices
#     # TODO: Figure out how to match orientations

#     return Tetrahedron(a, b, c, d), Tetrahedron(b, c, d, e)
# end
function flip14(tet::DelaunayTet, p::Int)
    a, b, c, d = tet.vertices

    t1 = DelaunayTet(a, b, c, p)
    t2 = DelaunayTet(a, b, d, p)
    t3 = DelaunayTet(a, c, d, p)
    t4 = DelaunayTet(b, d, c, p)

    if isassigned(tet.neighbors, 1)
        replace_neighbor!(tet.neighbors[1], tet, t1)
        t1.neighbors[1] = tet.neighbors[1]
    end
    if isassigned(tet.neighbors, 2)
        replace_neighbor!(tet.neighbors[2], tet, t2)
        t2.neighbors[2] = tet.neighbors[2]
    end
    if isassigned(tet.neighbors, 3)
        replace_neighbor!(tet.neighbors[3], tet, t3)
        t3.neighbors[3] = tet.neighbors[3]
    end
    if isassigned(tet.neighbors, 4)
        replace_neighbor!(tet.neighbors[4], tet, t4)
        t4.neighbors[4] = tet.neighbors[4]
    end

    t1.neighbors[2] = t2
    t1.neighbors[3] = t3
    t1.neighbors[4] = t4

    t2.neighbors[1] = t1
    t2.neighbors[3] = t3
    t2.neighbors[4] = t4

    t3.neighbors[1] = t1
    t3.neighbors[2] = t2
    t3.neighbors[4] = t4

    t4.neighbors[1] = t1
    t4.neighbors[2] = t2
    t4.neighbors[3] = t3

    return t1
end

# function flip41(tet1::Tetrahedron, tet2::Tetrahedron, tet3::Tetrahedron, tet4::Tetrahedron)
#     a, b, c, e = tet1.vertices
#     a, d, c, e = tet2.vertices
#     a, b, d, e = tet3.vertices
#     b, d, c, e = tet4.vertices
#     # TODO: Figure out how to match orientations

#     return Tetrahedron(a, b, c, d)
# end

# function flip44(tet1::Tetrahedron, tet2::Tetrahedron, tet3::Tetrahedron, tet4::Tetrahedron)
#     a, b, c, d = tet1.vertices
#     a, c, d, e = tet2.vertices
#     b, c, d, f = tet3.vertices
#     d, c, e, f = tet4.vertices
#     # TODO: Figure out how to match orientations

#     return Tetrahedron(a, b, c, e), Tetrahedron(a, b, e, d), Tetrahedron(b, c, e, f), Tetrahedron(b, e, d, f)
# end
