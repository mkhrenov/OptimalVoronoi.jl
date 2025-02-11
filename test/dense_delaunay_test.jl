using WeightedCVT
using SparseArrays
using GLMakie

d = 4
points = hcat([25 0 0 0; 0 25 0 0; 0 0 25 0], 2 .* rand(3, 20))
t = WeightedCVT.delaunay_tet(points)

@assert WeightedCVT.is_delaunay(t, points)

SMT = SparseMatrixCSC
r, dense_delaunay = WeightedCVT.condense_delaunay(t, points, SMT)
dense_voronoi = WeightedCVT.dual_complex(dense_delaunay)

# fig, ax, plot = WeightedCVT.viz(t, points)
# scatter!(dense_voronoi.vertices)
# scatter!(dense_voronoi.cell_centers)
# mesh!(Sphere(Point3(dense_voronoi.vertices[:, 1]), r[1]))

WeightedCVT.viz(dense_delaunay)
WeightedCVT.viz!(dense_voronoi, edge_color=:red)

boundary = WeightedCVT.get_boundary(dense_delaunay)
WeightedCVT.viz(boundary)