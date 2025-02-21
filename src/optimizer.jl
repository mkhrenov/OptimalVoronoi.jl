
# BFGS optimizer for voronoi tessellations, with line-search to ensure strict feasibility
function optimize_voronoi(f, g!, c!, Ω::F, x₀, Nx, Nc, μ; tol=1e-3, max_iters=1000) where {F}
    @assert Nx == length(x₀)
    α = 1.0

    xₖ₊₁ = copy(x₀)
    gₖ₊₁ = zeros(Nx)
    xₖ = copy(x₀)
    gₖ = zeros(Nx)
    xs = copy(x₀)

    c = zeros(Nc)           # Constraint evaluation
    y = copy(gₖ)            # Gradient secant
    p = copy(gₖ)            # Step direction
    s = copy(gₖ)            # Final step

    H = diagm(ones(Nx))     # Approximate Hessian inverse
    Hb = copy(H)
    outer_prod = zeros(Nx, Nx)

    points = reshape(x₀, 3, :)
    voronoi = bounded_voronoi(points, Ω)
    g!(gₖ₊₁, x₀, voronoi)
    for k in 1:max_iters
        # Get gradient at xₖ
        gₖ .= gₖ₊₁
        xₖ .= xₖ₊₁

        if maximum(abs, gₖ) < tol
            break
        end

        # Solve for direction pₖ
        mul!(p, H, gₖ)
        p .*= -1.0

        # Perform line-search
        α = linesearch(α, f, c!, c, Ω, p, xₖ, xs)

        # Update iterate
        @. s = α * p
        @. xₖ₊₁ = xₖ + s

        # Get gradient at xₖ₊₁, compute secant
        new_points = reshape(xₖ₊₁, 3, :)
        voronoi = bounded_voronoi(new_points, Ω)
        g!(gₖ₊₁, xₖ₊₁, voronoi)
        @. y = gₖ₊₁ - gₖ

        curvature = dot(s, y)
        if curvature ≤ 0.0
            y .+= (-curvature + 0.01) / dot(s, s) * s
            curvature = dot(s, y)
        end

        # Update inverse Hessian estimate
        # Hₖ .+= (curvature + dot(y, Hₖ, y)) * (s * s') / (curvature)^2 .- (Hₖ * y * s' + s * y' * Hₖ) / (curvature)
        Hb .= H
        mul!(outer_prod, s, s')
        H .+= ((curvature + dot(y, H, y)) / (curvature)^2) .* outer_prod
        mul!(outer_prod, y, s')
        mul!(H, Hb, outer_prod, -1.0 / curvature, 1.0)
        mul!(H, outer_prod', Hb, -1.0 / curvature, 1.0)


        α = clamp(α * 1.2, 0.0, 1.0)
        # if mod(k, 10) == 0
        #     println("$(maximum(abs, gₖ₊₁))")#$(round(f(xₖ₊₁), digits=5)), 
        # end
    end

    return xₖ
end

function linesearch(α, f, c!, c, Ω::F, p, x0, xs) where {F}
    for k in 1:32
        @. xs = x0 + α * p
        c!(c, xs, Ω)

        if minimum(c) > 0.1
            # if f(x) - f(x0) < 0.1
            return α
            # end
        end

        α *= sqrt(0.1)
    end
    println("Line-search failed, α: $α")

    return α
end