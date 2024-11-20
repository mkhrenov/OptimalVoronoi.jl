using Delaunay
using WeightedCVT
using BenchmarkTools
using Profile
using StaticArrays
using GeometryBasics
using GLMakie


points = rand(20, WeightedCVT.DIM)

tri = delaunay(points)
voronoi = WeightedCVT.delaunay_to_voronoi(tri)

tri_point_pairs = [
    (tri.simplices[i,j], tri.simplices[i,k])
    for i in 1:size(tri.simplices, 1)
    for j in 1:size(tri.simplices, 2)
    for k in (j+1):size(tri.simplices, 2)
]

tetras = [GeometryBasics.TetrahedronFace(tri.simplices[i, :]...) for i in 1:size(tri.simplices, 1)]
points = Makie.to_vertices(tri.points) # Use Makie to convert to Vector{Point3f}
m = GeometryBasics.Mesh(points, tetras) # create tetrahedra mesh
# Triangulate it, since Makie's mesh conversion currently doesn't handle tetrahedras itself 

# scatter(tri.points[:, 1], tri.points[:, 2], tri.points[:, 3], label = "Delaunay Points")
# plot([(tri.points[i,:], tri.points[j,:]) for (i,j) in tri_point_pairs])
# scatter!(voronoi.vertices[:, 1], voronoi.vertices[:, 2], voronoi.vertices[:, 3], label = "Voronoi Cell Vertices")
