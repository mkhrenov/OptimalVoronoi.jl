##### Flips #####

function to_vertices(tet1::DelaunayTet, tet2::DelaunayTet, points)
    n1 = get_neighbor_number(tet1, tet2)
    n2 = get_neighbor_number(tet2, tet1)
    b, c, d, a = get_commonface_and_opposite(tet1, n1)
    e = tet2.vertices[5-n2]

    return a, b, c, d, e
end

function flip23(tet1::DelaunayTet, tet2::DelaunayTet, points)
    a, b, c, d, e = to_vertices(tet1, tet2, points)

    @assert simplex_is(tet1, a, b, c, d)
    @assert simplex_is(tet2, b, c, d, e)

    t1 = DelaunayTet(a, b, c, e)
    @assert orient(a, b, c, e, points) > 0.0
    t2 = DelaunayTet(a, b, e, d)
    @assert orient(a, b, e, d, points) > 0.0
    t3 = DelaunayTet(a, e, c, d)
    @assert orient(a, e, c, d, points) > 0.0

    new_tets = SVector{3,DelaunayTet}(t1, t2, t3)
    ring_verts = SVector{3,Int}(d, c, b)

    for i in eachindex(ring_verts)
        t = get_neighbor_opposite(tet1, ring_verts[i])
        if t != tet1
            replace_neighbor!(t, tet1, new_tets[i])
            new_tets[i].neighbors[i] = t
        end

        t = get_neighbor_opposite(tet2, ring_verts[i])
        if t != tet2
            replace_neighbor!(t, tet2, new_tets[i])
            new_tets[i].neighbors[4] = t
        end
    end

    # Assign neighbors between new tetrahedra
    for i in eachindex(new_tets)
        for j in eachindex(new_tets)
            if i != j
                new_tets[i].neighbors[j] = new_tets[j]
            end
        end
    end

    return t1, t2, t3
end


function to_vertices(tet1::DelaunayTet, tet2::DelaunayTet, tet3::DelaunayTet, points)
    _, _, _, b = get_commonface_and_opposite(tet1, tet3)
    _, _, _, c = get_commonface_and_opposite(tet3, tet2)
    _, _, _, d = get_commonface_and_opposite(tet3, tet1)

    a = 0
    e = 0

    for v in tet1.vertices
        if (v != b) && (v != c) && (v != e)
            if orient(c, b, d, v, points) > 0
                a = v
            else
                e = v
            end
        end
    end

    return a, b, c, d, e
end

function flip32(tet1::DelaunayTet, tet2::DelaunayTet, tet3::DelaunayTet, points)
    a, b, c, d, e = to_vertices(tet1, tet2, tet3, points)

    t1 = DelaunayTet(a, b, c, d)
    @assert orient(a, b, c, d, points) > 0.0
    t2 = DelaunayTet(e, c, b, d)
    @assert orient(e, c, b, d, points) > 0.0

    @assert simplex_is(tet1, a, b, c, e)
    @assert simplex_is(tet2, a, b, e, d)
    @assert simplex_is(tet3, a, e, c, d)

    t1.neighbors[4] = t2
    t2.neighbors[4] = t1

    old_tets = SVector{3,DelaunayTet}(tet1, tet2, tet3)
    for i in eachindex(old_tets)
        t = get_neighbor_opposite(old_tets[i], a)
        if t != old_tets[i]
            replace_neighbor!(t, old_tets[i], t2)
            t2.neighbors[i == 1 ? 1 : (i == 2 ? 3 : 2)] = t
        end

        t = get_neighbor_opposite(old_tets[i], e)
        if t != old_tets[i]
            replace_neighbor!(t, old_tets[i], t1)
            t1.neighbors[i] = t
        end
    end

    return t1, t2
end

function flip14(tet::DelaunayTet, p::Int, points)
    a, b, c, d = tet.vertices

    t1 = DelaunayTet(a, b, c, p)
    @assert orient(a, b, c, p, points) > 0.0
    t2 = DelaunayTet(a, b, p, d)
    @assert orient(a, b, p, d, points) > 0.0
    t3 = DelaunayTet(a, p, c, d)
    @assert orient(a, p, c, d, points) > 0.0
    t4 = DelaunayTet(p, b, c, d)
    @assert orient(p, b, c, d, points) > 0.0

    new_tets = SVector{4,DelaunayTet}(t1, t2, t3, t4)

    for i in eachindex(new_tets)
        # Assign neighbors to neighbors of old tetrahedron
        if tet.neighbors[i] != tet
            replace_neighbor!(tet.neighbors[i], tet, new_tets[i])
            new_tets[i].neighbors[i] = tet.neighbors[i]
        end

        # Assign neighbors between new tetrahedra
        for j in eachindex(new_tets)
            if i != j
                new_tets[i].neighbors[j] = new_tets[j]
            end
        end
    end

    return t1, t2, t3, t4
end

# function flip41(tet1::Tetrahedron, tet2::Tetrahedron, tet3::Tetrahedron, tet4::Tetrahedron, points)
#     a, b, c, e = tet1.vertices
#     a, d, c, e = tet2.vertices
#     a, b, d, e = tet3.vertices
#     b, d, c, e = tet4.vertices
#     # TODO: Figure out how to match orientations

#     return Tetrahedron(a, b, c, d)
# end

# function flip44(tet1::Tetrahedron, tet2::Tetrahedron, tet3::Tetrahedron, tet4::Tetrahedron, points)
#     a, b, c, d = tet1.vertices
#     a, c, d, e = tet2.vertices
#     b, c, d, f = tet3.vertices
#     d, c, e, f = tet4.vertices
#     # TODO: Figure out how to match orientations

#     return Tetrahedron(a, b, c, e), Tetrahedron(a, b, e, d), Tetrahedron(b, c, e, f), Tetrahedron(b, e, d, f)
# end
