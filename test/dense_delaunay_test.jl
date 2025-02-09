using WeightedCVT
using SparseArrays
using GLMakie

d = 4
points = hcat([5 0 0 0; 0 5 0 0; 0 0 5 0], rand(3, 100))
t = WeightedCVT.delaunay_tet(points)

@assert WeightedCVT.is_delaunay(t, points)

SMT = SparseMatrixCSC
dense_delaunay = WeightedCVT.condense_delaunay(t, points, SMT)
dense_voronoi = WeightedCVT.dual_complex(dense_delaunay)

fig, ax, plot = WeightedCVT.viz(t, points)
scatter!(dense_voronoi.vertices)
scatter!(dense_voronoi.cell_centers)