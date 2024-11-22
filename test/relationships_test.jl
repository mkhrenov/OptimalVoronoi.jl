using WeightedCVT
using GLMakie
using BenchmarkTools
using Profile

const T0 = 0.1#273.15
const V = 5.0 # mm / s
const P = 1000.0 # W
const k = 45.e-3 # W / mm K 
const ρ = 7850.0e-6 # g / mm³
const cₚ = 420.0e-3 # J / g K
const α = k / ρ / cₚ
const nx = 100
const ny = 100
const nz = 100

T(x) = T0 + (P / (2π * k * √((x[1] - nx ÷ 2)^2 + (x[2] - ny ÷ 2)^2 + x[3]^2))) * exp(-V * (√((x[1] - nx ÷ 2)^2 + (x[2] - ny ÷ 2)^2 + x[3]^2) - (x[1] - nx ÷ 2)) / (2α))
N = 300

domain = ones(Int, nx, ny, nz)
for index in CartesianIndices(domain)
    dist_to_c = (index[1] - nx ÷ 2)^2 + (index[2] - ny ÷ 2)^2 + (index[3])^2
    if dist_to_c > 60^2
        domain[index] = 0
    end
end
min_dist = zeros(size(domain))

init_points = Float64.(rand(1:ny, 3, N))
points = WeightedCVT.centroidal_voronoi(domain, init_points, T)

volumes, adjacency, separations, interfaces, air_area, base_area, meshes = WeightedCVT.structure(domain, points)

meshes = [mesh for mesh in meshes if !isempty(mesh)]
fig = GLMakie.mesh(meshes[1])
GLMakie.mesh!.(meshes)
scatter!(points, color=:red)
display(fig)