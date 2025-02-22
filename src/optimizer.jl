
# BFGS optimizer for voronoi tessellations, with line-search to ensure strict feasibility
function optimize_voronoi(f, g!, c!, Ω::F, x₀, Nx, Nc;
    tol = 1e-3, max_iters = 1000, α = 1.0, β₁ = 0.9, β₂ = 0.999, ϵ = 1e-8) where {F}
    @assert Nx == length(x₀)

    α_init = α
    x = copy(x₀)
    g = zeros(Nx)
    m = zeros(Nx)
    v = zeros(Nx)
    dx = copy(x)

    p = zeros(Nx)

    xs = copy(x₀)
    c = zeros(Nc)           # Constraint evaluation

    for k in 1:max_iters
        points = reshape(x, 3, :)
        voronoi = bounded_voronoi(points, Ω)

        # Compute gradient
        g!(g, x, voronoi)
        penalty!(g, c!, c, Ω, x, 1e-6)

        # @show g

        # Update biased moment estimates
        @. m = β₁ * m + (1 - β₁) * g
        @. v = β₂ * v + (1 - β₂) * g^2

        # Compute step direction
        @. p = -1 * m / (1 - β₁^k) / (√(v / (1 - β₂^k)) + ϵ)

        # Perform line-search to ensure feasibility
        α = linesearch(α, f, c!, c, Ω, p, x, xs)


        # Update iterate
        # dx .= x
        @. x += α * p
        # dx .-= x

        project_to_volume!(x, Ω)

        if maximum(abs, g) ≤ tol
            break
        end

        α = clamp(α * 1.2, 0.0, α_init)
        println("$(maximum(abs, g))")#$(round(f(xₖ₊₁), digits=5)), 
    end

    return x
end

function project_to_volume!(x, Ω::F) where {F}
    points = reshape(x, 3, :)
    n_points = size(points, 2)
    g = zeros(3)

    for i in 1:n_points
        p_i = view(points, :, i)

        while Ω(p_i) > -0.5
            ForwardDiff.gradient!(g, Ω, p_i)
            p_i .-= (Ω(p_i) + 0.51) .* g
        end
    end
end

function penalty!(g, c!, c, Ω::F, x, μ) where {F}
    Nc = length(c)
    points = reshape(x, 3, :)
    gm = reshape(g, 3, :)
    c!(c, x, Ω)

    for i in 1:Nc
        gm[:, i] .-= μ / c[i] * ForwardDiff.gradient(Ω, points[:, i])
    end
end

function linesearch(α, f, c!, c, Ω::F, p, x0, xs) where {F}
    while α ≥ 1e-16
        @. xs = x0 + α * p
        c!(c, xs, Ω)

        if maximum(c) < 0.0
            return α
        end
        println("Overshoot of $(minimum(c)), α=$α")

        α *= sqrt(0.1)
    end
    println("Line-search failed, α: $α")

    return α
end