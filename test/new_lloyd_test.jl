using WeightedCVT
using SparseArrays
using GLMakie
using Profile
using Random

sdf = zeros(30, 40, 30)
for cidx in CartesianIndices(sdf)
    if cidx[3] >= 10
        sdf[cidx] = min(
            √((cidx[1] - 10)^2 + (cidx[2] - 20)^2 + (cidx[3] - 10)^2) - 7,
            √((cidx[1] - 20)^2 + (cidx[2] - 20)^2 + (cidx[3] - 10)^2) - 7,
            √((cidx[1] - 15)^2 + (cidx[2] - 20)^2 + (cidx[3] - 10)^2) - 7,
        )
    else
        sdf[cidx] = min(
            max(abs(cidx[1] - 10), abs(cidx[2] - 20), abs(cidx[3] - 10)) - 9,
            max(abs(cidx[1] - 20), abs(cidx[2] - 20), abs(cidx[3] - 10)) - 9,
            max(abs(cidx[1] - 15), abs(cidx[2] - 20), abs(cidx[3] - 10)) - 9,
        )
    end
end
Ω = WeightedCVT.Ω_from_array(sdf)

points = zeros(3, 100)
x = zeros(3)
for k in 1:100
    x[1] = 30 * rand()
    x[2] = 40 * rand()
    x[3] = 30 * rand()

    while Ω(x) ≥ -1
        x[1] = 30 * rand()
        x[2] = 40 * rand()
        x[3] = 30 * rand()
    end

    points[:, k] .= x
end
pcopy = copy(points)
scatter(points)

points .= pcopy
@time voronoi = WeightedCVT.lloyd(points, Ω);

points .= pcopy
@time voronoi = WeightedCVT.lloyd(points, Ω);

# points .= pcopy
# @profview voronoi = WeightedCVT.lloyd(points, Ω);

# points .= pcopy
# @profview_allocs voronoi = WeightedCVT.lloyd(points, Ω) sample_rate=0.0001

voronoi.vertices .-= 0.6
voronoi.cell_centers .-= 0.6
WeightedCVT.viz(voronoi)
volume!(sdf, algorithm=:iso, isovalue=0, isorange=0.1, alpha=0.1)

# bd = WeightedCVT.get_boundary(voronoi)
# WeightedCVT.viz(bd)