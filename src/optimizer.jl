
# BFGS optimizer for voronoi tessellations, with line-search to ensure strict feasibility
function optimize_voronoi(f, g!, c!, Ω::F, x₀, Nx, Nc;
    tol=1e-3, max_iters=1000, α=1.0, β₁=0.9, β₂=0.999, ϵ=1e-8) where {F}
    @assert Nx == length(x₀)

    α_init = α

    xₖ₊₁ = copy(x₀)
    gₖ₊₁ = zeros(Nx)
    xₖ = copy(x₀)
    gₖ = zeros(Nx)
    xs = copy(x₀)

    c = zeros(Nc)           # Constraint evaluation
    y = copy(gₖ)            # Gradient secant
    p = copy(gₖ)            # Step direction
    s = copy(gₖ)            # Final step
    r = copy(y)

    gb = copy(gₖ)

    B = diagm(ones(Nx))      # Approximate Hessian inverse
    Bk = copy(B)
    outer_prod1 = zeros(Nx, Nx)
    outer_prod2 = zeros(Nx, Nx)

    # Evaluate first gradient
    points = reshape(xₖ₊₁, 3, :)
    @show J0 = f(xₖ₊₁, Ω)
    g!(gₖ₊₁, xₖ₊₁, Ω)
    penalty_gradient!(gₖ₊₁, c!, c, Ω, xₖ₊₁, 0.1)
    ρ = 1.0

    for b in 3.0:0.5:6.0
        μ = (0.1)^b
        println("######### $b ##########")

        B .= 0.0
        B[diagind(B)] .= 1.0
        Bk .= B

        # Get gradient at xₖ₊₁
        g!(gₖ₊₁, xₖ₊₁, Ω)
        penalty_gradient!(gₖ₊₁, c!, c, Ω, xₖ₊₁, μ)

        for k in 1:max_iters
            # Get gradient at xₖ
            gₖ .= gₖ₊₁
            xₖ .= xₖ₊₁

            if maximum(abs, gₖ) < tol
                break
            end

            # Solve for direction pₖ
            B .+= B'
            B ./= 2.0
            Bk .= B
            B_factorized = LinearAlgebra.cholesky!(Hermitian(B), check=false)

            if !issuccess(B_factorized)
                B .= Bk
                B .= 0.0
                B[diagind(B)] .+= ρ
                # @show ρ
                ρ *= 2.0
                # α = α_init
                continue
                # else
                #     ρ = 1.0
            end

            ldiv!(p, B_factorized, gₖ)
            p .*= -1.0
            B .= Bk

            # Perform line-search to ensure feasibility
            α = linesearch(α, f, g!, c!, c, Ω, p, xₖ, xs, gₖ, gb, μ, J0)

            if α < 0.0
                B .= 0.0
                B[diagind(B)] .+= ρ
                @show ρ
                ρ *= 2.0
                α = α_init
                continue
            else
                ρ = 1.0
            end

            # @show α
            # Update iterate
            @. s = α * p
            @. xₖ₊₁ = xₖ + s

            # Get gradient at xₖ₊₁, compute secant
            g!(gₖ₊₁, xₖ₊₁, Ω)
            penalty_gradient!(gₖ₊₁, c!, c, Ω, xₖ₊₁, μ)
            @. y = gₖ₊₁ - gₖ

            sTy = dot(s, y)
            sBs = dot(s, B, s)
            # @show sBs
            # @show sTy

            # Compute curvature damping
            θ = (sTy ≥ 0.2 * sBs) ? 1.0 : (0.8 * sBs) / (sBs - sTy)
            # @show θ

            # Update damped secant
            # r = θ * y + (1 - θ) * B * s
            mul!(r, B, s)
            r .*= (1.0 - θ)
            r .+= θ .* y

            sTr = dot(s, r)

            # Update Hessian estimate
            # B = B - (B * s * s' * B) / (s' * B * s) + (r * r') / (s' * r)
            Bk .= B
            mul!(outer_prod1, s, s')
            mul!(outer_prod2, Bk, outer_prod1)
            mul!(B, outer_prod2, Bk, -1.0 / sBs, 1.0)
            mul!(B, r, r', 1.0 / sTr, 1.0)

            α = clamp(α * 1.2, 0.0, α_init)
            println("$(maximum(abs, gₖ₊₁))")#$(round(f(xₖ₊₁), digits=5)), 
        end

    end

    return xₖ
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

function penalty(c!, c, Ω::F, x, μ) where {F}
    Nc = length(c)
    J = 0.0
    c!(c, x, Ω)

    for i in 1:Nc
        J -= μ * log(-(c[i] + 0.5))
    end
    return J
end

function penalty_gradient!(g, c!, c, Ω::F, x, μ) where {F}
    Nc = length(c)
    points = reshape(x, 3, :)
    gm = reshape(g, 3, :)
    c!(c, x, Ω)

    for i in 1:Nc
        gm[:, i] .-= (μ / (c[i] + 0.5)) * ForwardDiff.gradient(Ω, points[:, i])
    end
end

function linesearch(α, f, g!, c!, c, Ω::F, p, x0, xs, g, gb, μ, J0) where {F}
    while α ≥ 1e-8
        @. xs = x0 + α * p
        c!(c, xs, Ω)

        if (maximum(c) < -0.51)
            if (f(xs, Ω) + penalty(c!, c, Ω, xs, μ) ≤ (f(x0, Ω) + penalty(c!, c, Ω, x0, μ)) * 1.1)# - dot(p, g))  # (1e-4) * (1/α) * 

            g!(gb, xs, Ω)
            penalty_gradient!(gb, c!, c, Ω, xs, μ)

            # if -dot(p, gb) ≤ -0.9 * dot(p, g)
            # if dot(gb, gb) ≤ dot(g, g)
                return α
            # end
            # else
            #     println("Wolfe II violated")
            # end
            # else
            #     println("Wolfe I violated")
            #     # @show f(xs, Ω) + penalty(c!, c, Ω, xs, μ)
            #     # @show f(x0, Ω) + penalty(c!, c, Ω, x0, μ) + (1e-4) * α * dot(p, g)
            #     # throw(ErrorException("Wolfe!"))
            end
        else
            println("Overshoot of $(maximum(c)), α=$α")
        end

        α *= sqrt(0.1)
    end
    println("Line-search failed, α: $α")

    return -1.0
end