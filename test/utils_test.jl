using WeightedCVT
using StaticArrays
using BenchmarkTools
# using Profiling

N = WeightedCVT.DIM+1

points = SMatrix{N,D}(rand(N,D))
c = WeightedCVT.circumcenter(points)