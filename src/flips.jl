
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

# function flip14(tet1::Tetrahedron, p::Point)
#     a, b, c, d = tet1.vertices

#     return Tetrahedron(a, b, c, p), Tetrahedron(a, d, c, p), Tetrahedron(a, b, d, p), Tetrahedron(b, d, c, p)
# end

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
