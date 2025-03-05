using OptimalVoronoi

points = hcat([5 0 0 0; 0 5 0 0; 0 0 5 0], rand(3, 100))
t = delaunay_tet(points)

@assert is_delaunay(t, points)

viz(t, points)