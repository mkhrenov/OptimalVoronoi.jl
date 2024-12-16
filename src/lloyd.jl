
function centroidal_voronoi(domain, initial_points, f::F; max_iters=50, over_relax=1.5, tol=0.1) where {F}
    old_points = copy(initial_points)
    new_points = copy(initial_points)
    N = size(initial_points, 2)

    for i in 1:max_iters
        voronoi!(domain, old_points)
        get_centroids!(domain, new_points, f)

        @. new_points = (new_points - old_points) * over_relax + old_points
        if mapreduce((x, y) -> (x - y)^2, +, new_points, old_points) / N < tol
            println(i)
            break
        end
        old_points .= new_points
    end

    return new_points
end

function centroidal_voronoi(domain, initial_points; max_iters=50, over_relax=1.5, tol=0.1)
    old_points = copy(initial_points)
    new_points = copy(initial_points)
    N = size(initial_points, 2)

    for i in 1:max_iters
        voronoi!(domain, old_points)
        get_centroids!(domain, new_points)

        @. new_points = (new_points - old_points) * over_relax + old_points
        if mapreduce((x, y) -> (x - y)^2, +, new_points, old_points) / N < tol
            println(i)
            break
        end
        old_points .= new_points
    end

    return new_points
end
