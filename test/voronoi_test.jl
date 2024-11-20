using WeightedCVT
using GLMakie
using BenchmarkTools

nx = ny = nz = 100

domain = ones(Int, nx, ny, nz)
for index in CartesianIndices(domain)
    dist_to_c = (index[1] - nx ÷ 2)^2 + (index[2] - ny ÷ 2)^2 + (index[3] - nz ÷ 2)^2
    if dist_to_c > 50^2
        domain[index] = 0
    end
end
min_dist = zeros(size(domain))

points = rand(1:ny, 3, 10)

WeightedCVT.voronoi!(domain, min_dist, points)

fig = volume(domain)#, aspect_ratio = :equal)
scatter!(points, label = nothing)
display(fig)
