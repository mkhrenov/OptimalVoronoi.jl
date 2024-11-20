
function weighed_cvt(points, boundary, f::F) where {F}
 
    while true
        triangulation = delaunay(points)
        tesselation = delaunay_to_voronoi(triangulation, boundary)
        points = centroids(tesselation)
        # TODO over-relaxation
    end
end

