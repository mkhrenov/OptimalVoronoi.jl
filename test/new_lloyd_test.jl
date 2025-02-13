using WeightedCVT
using SparseArrays
using GLMakie
using Profile

sdf = zeros(20, 20, 20)
for cidx in CartesianIndices(sdf)
    sdf[cidx] = √((cidx[1] - 10)^2 + (cidx[2] - 10)^2 + (cidx[3] - 10)^2) - 8
end
Ω = WeightedCVT.Ω_from_array(sdf)

points = 8 .* rand(3, 25) .+ 6
pcopy = copy(points)
@time voronoi = WeightedCVT.lloyd(points, Ω);

points .= pcopy
@time voronoi = WeightedCVT.lloyd(points, Ω);

# points .= pcopy
# @profview voronoi = WeightedCVT.lloyd(points, Ω);

voronoi.vertices .-= 0.5
WeightedCVT.viz(voronoi)
volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.05)

# bd = WeightedCVT.get_boundary(voronoi)
# WeightedCVT.viz(bd)