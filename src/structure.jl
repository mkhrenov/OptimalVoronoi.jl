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

    xmin, xmax, ymin, ymax, zmin, zmax = cell_bounds(domain)

    neighbor_incs::Vector{CartesianIndex} = [CartesianIndex((i == d ? inc : 0 for i in 1:D)...) for d in 1:D for inc in (1, -1)]

    # Compute volumes and identify neighbors
    for index in CartesianIndices(domain)
        cell = domain[index]
        if cell == 0
            continue
        end

        volumes[cell] += 1.0

        # Check for neighbors
        for neighbor_inc in neighbor_incs
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
    algorithm = MarchingCubes(iso=0.5)
    mask = BitArray(undef, dims)
    point_vec = Vector{Vector{Point3{Float64}}}(undef, N)
    trian_vec = Vector{Vector{TriangleFace{Int}}}(undef, N)
    meshes = Vector{GeometryBasics.Mesh}(undef, N)
    for cell in 1:N
        mask .= (domain .== cell)

        vertices, triangles = isosurface(
            view(mask, xmin[cell]:xmax[cell], ymin[cell]:ymax[cell], zmin[cell]:zmax[cell]),
            algorithm, (xmin[cell]:xmax[cell]) .* 1.0, (ymin[cell]:ymax[cell]) .* 1.0, (zmin[cell]:zmax[cell]) .* 1.0
        )

        vertex_points = [Point3(v[1], v[2], v[3]) for v in vertices]
        triangle_faces = [TriangleFace(t...) for t in triangles]
        point_vec[cell] = vertex_points
        trian_vec[cell] = triangle_faces
        air_area[cell] = area(vertex_points, triangle_faces)

        m = GeometryBasics.Mesh(vertex_points, triangle_faces)
        meshes[cell] = m
    end

    # Compute areas for each interface
    for (i, j, v) in zip(findnz(adjacency)...)
        interfaces[i, j] = masked_area(point_vec[i], trian_vec[i], domain, j)
        air_area[i] -= interfaces[i, j]
    end

    # Take midpoint between area estimates
    interfaces .+= transpose(interfaces)
    interfaces ./= 2.0

    separations = sparse(separations)
    interfaces = sparse(interfaces)

    return volumes, adjacency, separations, interfaces, air_area, base_area, meshes
end

function masked_area(vertices::Vector{Point3{Float64}}, triangles::Vector{TriangleFace{Int}}, domain, n)
    Nx, Ny, Nz = size(domain)
    cum_area = 0.0

    for triangle in triangles
        in_interface = true

        for p in triangle
            point = vertices[p]
            i, j, k = round(Int, point[1]), round(Int, point[2]), round(Int, point[3])

            if (
                (domain[min(i + 1, Nx), max(j - 1, 1), max(k - 1, 1)] == n) || (domain[min(i + 1, Nx), max(j - 1, 1), k] == n) || (domain[min(i + 1, Nx), max(j - 1, 1), min(k + 1, Nz)] == n) ||
                (domain[min(i + 1, Nx), j, max(k - 1, 1)] == n) || (domain[min(i + 1, Nx), j, k] == n) || (domain[min(i + 1, Nx), j, min(k + 1, Nz)] == n) ||
                (domain[min(i + 1, Nx), min(j + 1, Ny), max(k - 1, 1)] == n) || (domain[min(i + 1, Nx), min(j + 1, Ny), k] == n) || (domain[min(i + 1, Nx), min(j + 1, Ny), min(k + 1, Nz)] == n) ||
                (domain[i, max(j - 1, 1), max(k - 1, 1)] == n) || (domain[i, max(j - 1, 1), k] == n) || (domain[i, max(j - 1, 1), min(k + 1, Nz)] == n) ||
                (domain[i, j, max(k - 1, 1)] == n) || (domain[i, j, k] == n) || (domain[i, j, min(k + 1, Nz)] == n) ||
                (domain[i, min(j + 1, Ny), max(k - 1, 1)] == n) || (domain[i, min(j + 1, Ny), k] == n) || (domain[i, min(j + 1, Ny), min(k + 1, Nz)] == n) ||
                (domain[max(i - 1, 1), max(j - 1, 1), max(k - 1, 1)] == n) || (domain[max(i - 1, 1), max(j - 1, 1), k] == n) || (domain[max(i - 1, 1), max(j - 1, 1), min(k + 1, Nz)] == n) ||
                (domain[max(i - 1, 1), j, max(k - 1, 1)] == n) || (domain[max(i - 1, 1), j, k] == n) || (domain[max(i - 1, 1), j, min(k + 1, Nz)] == n) ||
                (domain[max(i - 1, 1), min(j + 1, Ny), max(k - 1, 1)] == n) || (domain[max(i - 1, 1), min(j + 1, Ny), k] == n) || (domain[max(i - 1, 1), min(j + 1, Ny), min(k + 1, Nz)] == n)
            )
                continue
            else
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

function cell_bounds(domain)
    Nx, Ny, Nz = size(domain)
    Nc = maximum(domain)

    xmin = Nx .* ones(Int, Nc)
    ymin = Ny .* ones(Int, Nc)
    zmin = Nz .* ones(Int, Nc)
    xmax = zeros(Int, Nc)
    ymax = zeros(Int, Nc)
    zmax = zeros(Int, Nc)

    for x in 1:Nx
        for y in 1:Ny
            for z in 1:Nz
                cell = domain[x, y, z]

                xmin[cell] = min(max(x - 1, 1), xmin[cell])
                ymin[cell] = min(max(y - 1, 1), ymin[cell])
                zmin[cell] = min(max(z - 1, 1), zmin[cell])

                xmax[cell] = max(min(x + 1, Nx), xmax[cell])
                ymax[cell] = max(min(y + 1, Ny), ymax[cell])
                zmax[cell] = max(min(z + 1, Nz), zmax[cell])
            end
        end
    end

    return xmin, xmax, ymin, ymax, zmin, zmax
end
