using WeightedCVT
using SparseArrays
using GLMakie

sdf = zeros(20, 20, 20)
for cidx in CartesianIndices(sdf)
    sdf[cidx] = √((cidx[1] - 10)^2 + (cidx[2] - 10)^2 + (cidx[3] - 10)^2) - 8
end

d = 4
points = hcat([100 0 0 0; 0 100 0 0; 0 0 100 0], 8 .* rand(3, 25) .+ 6)
@time t = WeightedCVT.delaunay_tet(points)

@assert WeightedCVT.is_delaunay(t, points)

SMT = SparseMatrixCSC
@time r, dense_delaunay = WeightedCVT.condense_delaunay(t, points, SMT)
@time dense_voronoi = WeightedCVT.dual_complex(dense_delaunay)

Ω = WeightedCVT.Ω_from_array(sdf)
to_drop = .!WeightedCVT.in_Ω(dense_voronoi.vertices, Ω)

np = WeightedCVT.bound_voronoi(dense_voronoi, dense_delaunay, Ω)
WeightedCVT.viz(np)
volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.5, alpha=0.05)