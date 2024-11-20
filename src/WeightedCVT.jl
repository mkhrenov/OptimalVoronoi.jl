module WeightedCVT

using Delaunay
# using GeometryBasics

const DIM = 3

include("utils.jl")
include("voronoi.jl")
include("centroids.jl")
include("lloyd.jl")

end # module WeightedCVT
