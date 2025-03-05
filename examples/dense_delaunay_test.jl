using OptimalVoronoi
using SparseArrays
using GLMakie

points = hcat([25 0 0 0; 0 25 0 0; 0 0 25 0], 2 .* rand(3, 20))
t = delaunay_tet(points)

@assert is_delaunay(t, points)

r, dense_delaunay = condense_delaunay(t, points)
dense_voronoi = dual_complex(dense_delaunay)

viz(dense_delaunay)
viz!(dense_voronoi, edge_color=:red)