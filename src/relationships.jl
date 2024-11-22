"""
    Calculate the mesh structure from a voxelized Voronoi tesselation,
    specified using voxel-grid `domain` and generator points `points`
"""
function structure(domain::AbstractArray{<:Int}, points::AbstractMatrix)
    dims = size(domain)
    D, N = size(points)
    volumes = zeros(N)
    adjacency = falses(N, N)
    separations = zeros(N, N)
    interfaces = zeros(N, N)
    air_area = zeros(N)
    base_area = zeros(N)

    neighbor_incs = [CartesianIndex((i == d ? inc : 0 for i in 1:D)...) for d in 1:D for inc in (1, -1)]

    # Compute volumes and identify neighbors
    for index in CartesianIndices(domain)
        cell = domain[index]
        if cell == 0
            continue
        end

        volumes[cell] += 1.0

        # Check for neighbors
        for neighbor_inc::CartesianIndex in neighbor_incs
            neighbor_index = index + neighbor_inc
            # Make sure this neighbor exists (we aren't out of bounds)
            if !checkbounds(Bool, domain, neighbor_index)
                base_area[cell] += 1.0
                continue
            end

            neighbor_cell::Int = domain[neighbor_index]
            if (neighbor_cell != 0) && (!adjacency[cell, neighbor_cell]) && (neighbor_cell != cell)
                adjacency[cell, neighbor_cell] = true
                adjacency[neighbor_cell, cell] = true
            end
        end
    end
    adjacency = sparse(adjacency)

    # Compute separations for neighbors
    dr = zeros(D)
    for (i, j, v) in zip(findnz(adjacency)...)
        @. dr = points[i] - points[j]
        separations[i, j] = norm(dr)
    end

    # Mesh each cell
    algorithm = MarchingCubes(iso = 0.5)
    mask = BitArray(undef, dims)
    point_vec = Vector{Vector{Point3{Float64}}}(undef, N)
    trian_vec = Vector{Vector{TriangleFace}}(undef, N)
    meshes = Vector{GeometryBasics.Mesh}(undef, N)
    for cell in 1:N
        map!(v -> (v == cell) ? true : false, mask, domain)

        # TODO : Optimize marching cubes (shrink domain)
        vertices, triangles = isosurface(mask, algorithm, 1:dims[1], 1:dims[2], 1:dims[3])

        vertex_points = [Point3(v...) for v in vertices]
        triangle_faces = [TriangleFace(t...) for t in triangles]
        point_vec[cell] = vertex_points
        trian_vec[cell] = triangle_faces
        air_area[cell] = area(vertex_points, triangle_faces)

        m = GeometryBasics.Mesh(vertex_points, triangle_faces)
        meshes[cell] = m
    end

    # Compute areas for each interface
    for (i, j, v) in zip(findnz(adjacency)...)
        mask .= false

        # Compute boundary layer for cell pair (i, j)
        for index in CartesianIndices(domain)
            if domain[index] != i
                continue
            end

            for neighbor_inc::CartesianIndex in neighbor_incs
                neighbor_index = index + neighbor_inc
                # Make sure this neighbor exists (we aren't out of bounds)
                if !checkbounds(Bool, domain, neighbor_index)
                    continue
                end

                if domain[neighbor_index] == j
                    mask[index] = true
                    mask[neighbor_index] = true
                end
            end
        end

        # Cull triangles that aren't part of the interface (i, j)
        vertex_points = point_vec[i]
        triangle_faces = trian_vec[i]
        interfaces[i, j] = masked_area(vertex_points, triangle_faces, mask)
        air_area[i] -= interfaces[i, j]
    end

    # Take midpoint between area estimates
    @. interfaces = (interfaces + interfaces') / 2.0

    separations = sparse(separations)
    interfaces = sparse(interfaces)

    return volumes, adjacency, separations, interfaces, air_area, base_area, meshes
end

function masked_area(vertices::AbstractVector{<:Point3}, triangles::AbstractVector{<:TriangleFace}, mask)
    cum_area = 0.0

    for triangle::TriangleFace in triangles
        in_interface = true
        for p in triangle
            point::Point3{Float64} = vertices[p]
            index = CartesianIndex(round.(Int, point)...)

            if !mask[index]
                in_interface = false
                break
            end
        end

        if in_interface
            cum_area += area(vertices, triangle)
        end
    end

    return cum_area
end
