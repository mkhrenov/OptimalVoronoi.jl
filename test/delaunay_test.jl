using WeightedCVT

d = 4
points = hcat([5 0 0 0; 0 5 0 0; 0 0 5 0], rand(3, 100))
t = WeightedCVT.delaunay_tet(points)

@assert WeightedCVT.is_delaunay(t, points)

WeightedCVT.viz(t, points)