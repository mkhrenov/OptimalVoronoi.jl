using WeightedCVT
using SparseArrays
using GLMakie

points = hcat([25 0 0 0; 0 25 0 0; 0 0 25 0], 2 .* rand(3, 20))
t = WeightedCVT.delaunay_tet(points)

@assert WeightedCVT.is_delaunay(t, points)

r, dense_delaunay = WeightedCVT.condense_delaunay(t, points)
dense_voronoi = WeightedCVT.dual_complex(dense_delaunay)

WeightedCVT.viz(dense_delaunay)
WeightedCVT.viz!(dense_voronoi, edge_color=:red)