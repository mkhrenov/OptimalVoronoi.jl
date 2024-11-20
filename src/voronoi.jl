
struct VoronoiTesselation{T}
    dim::Int
    generator_points::Matrix{T}
    vertices::Matrix{T}
    cell_vertices::Vector{Vector{Int}}
    cell_faces::Vector{Vector{Int}}
    face_vertices::Vector{Vector{Int}}
    face_areas::Vector{T}
    face_separations::Vector{T}
    # exterior????

    function VoronoiTesselation{T}(dim, n_cells, n_vertices, n_faces) where {T}
        generator_points = zeros(T, (n_cells, dim))
        vertices = zeros(T, (n_vertices, dim))
        cell_vertices = Vector{Vector{Int}}()
        cell_faces = Vector{Vector{Int}}()
        face_vertices = Vector{Vector{Int}}()
        face_areas = zeros(T, (n_faces,))
        face_separations = zeros(T, (n_faces,))

        return new{T}(dim,
            generator_points,
            vertices,
            cell_vertices,
            cell_faces,
            face_vertices,
            face_areas,
            face_separations)
    end
end

# For interior points, 
# each Delaunay point maps to a Voronoi cell
# each Delaunay edge maps to a Voronoi face
# each Delaunay simplex face maps to a Voronoi edge
# each Delaunay simplex volume maps to a Voronoi vertex

function delaunay_to_voronoi(DT::Triangulation) #, boundary::Mesh)
    simplices = DT.simplices
    generator_points = DT.points
    edge_connections = DT.neighbors
    face_connections = DT.vertex_neighbor_vertices
    n_cells, dim = size(generator_points)
    n_vertices = size(simplices, 1) + 0 # Add boundary
    n_faces = sum(face_connections) ÷ 2

    voronoi = VoronoiTesselation{Float64}(dim, n_cells, n_vertices, n_faces)
    voronoi.generator_points .= generator_points

    # Loop over simplices, find circumcenters to get Voronoi vertices in the interior
    for i in 1:size(simplices, 1)
        points = zeros(SMatrix{0, DIM})
        for j in 1:(DIM+1)
            vertex = simplices[i, j]
            new_point = SVector{DIM}(@view generator_points[vertex, :])
            points = vcat(points, transpose(new_point))
        end

        circum = circumcenter(points)
        voronoi.vertices[i, :] .= circum
    end

    

    # Connect voronoi vertices with edges (using Simplex neighbors)
    for edge in edge_connections
    end

    # Connect voronoi vertices to get faces (using Delaunay edges)
    for face in face_connections
    end

    # Use faces to get surface areas and Voronoi cell sub-simplices

    # Outer boundaries??????

    return voronoi
end
