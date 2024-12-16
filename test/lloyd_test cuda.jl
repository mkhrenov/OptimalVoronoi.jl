using WeightedCVT
using GLMakie
using BenchmarkTools
using Profile
using CUDA

nx = ny = nz = 100
N = 10

domain = ones(Int, nx, ny, nz)
for index in CartesianIndices(domain)
    dist_to_c = (index[1] - nx ÷ 2)^2 + (index[2] - ny ÷ 2)^2 + (index[3] - nz ÷ 2)^2
    if dist_to_c > 60^2
        domain[index] = 0
    end
end
points = Float64.(rand(1:ny, 3, N))

domain = cu(domain)
points = cu(points)

points = WeightedCVT.centroidal_voronoi(domain, points)

domain = Array(domain)
points = Array(points)

fig = volume(domain)
scatter!(points, color=:red)
display(fig)
