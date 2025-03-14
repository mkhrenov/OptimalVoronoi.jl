
mutable struct BFGSOptimizerState{F,G,C,ST,VT,MT}
    ρ::ST
    μ::ST
    α::ST
    lim::ST

    max_iters::Int
    verbosity::Int

    f_abs_tol::ST
    g_abs_tol::ST
    f_rel_tol::ST
    g_rel_tol::ST

    f::F
    g!::G
    Ω::C

    fₖ::ST
    fₖ₊₁::ST
    xₖ::VT
    xₖ₊₁::VT
    gₖ::VT
    gₖ₊₁::VT
    cₖ::VT
    cₖ₊₁::VT


    Bₖ::MT
    Bₒ::MT
    o_prod_1::MT
    o_prod_2::MT

    y::VT
    p::VT
    s::VT
    r::VT

    function BFGSOptimizerState(
        f::F, g!::G, Ω::C, Nx, Nc; μ₀=0.1, α₀=1.0, max_iters=10_000, verbosity=4,
        f_abs_tol=1e-6, g_abs_tol=1e-4, f_rel_tol=1e-6, g_rel_tol=1e-4,
        ST=Float64, VT=Vector{ST}, MT=Matrix{ST}) where {F,G,C}

        xₖ = zeros(Nx)
        xₖ₊₁ = zeros(Nx)
        gₖ = zeros(Nx)
        gₖ₊₁ = zeros(Nx)
        cₖ = zeros(Nc)
        cₖ₊₁ = zeros(Nc)

        Bₖ = zeros(Nx, Nx)
        Bₒ = zeros(Nx, Nx)
        o_prod_1 = zeros(Nx, Nx)
        o_prod_2 = zeros(Nx, Nx)

        y = zeros(Nx)
        p = zeros(Nx)
        s = zeros(Nx)
        r = zeros(Nx)

        return new{F,G,C,ST,VT,MT}(
            0.0, μ₀, α₀, 0.5,
            max_iters, verbosity,
            f_abs_tol, g_abs_tol, f_rel_tol, g_rel_tol,
            f, g!, Ω,
            0.0, 0.0,
            xₖ, xₖ₊₁, gₖ, gₖ₊₁, cₖ, cₖ₊₁,
            Bₖ, Bₒ, o_prod_1, o_prod_2,
            y, p, s, r
        )
    end
end

function initialize!(o::BFGSOptimizerState, x₀)
    o.xₖ .= x₀
    o.xₖ₊₁ .= x₀

    o.fₖ = opt_objective(o, o.xₖ, o.cₖ)
    o.fₖ₊₁ = o.fₖ

    opt_gradient!(o, o.gₖ, o.xₖ, o.cₖ)
    o.gₖ₊₁ .= o.gₖ

    reset_hessian!(o)
end

function reset_hessian!(o::BFGSOptimizerState)
    o.ρ = 0.0
    o.α = 1.0
    o.Bₖ .= 0.0
    o.Bₖ[diagind(o.Bₖ)] .= 1.0 * norm(o.gₖ)
    o.Bₒ .= o.Bₖ
end

function opt_objective(o::BFGSOptimizerState, x, c)
    # Compute main objective
    obj = o.f(x)

    # Add penalties SDF
    points = reshape(x, 3, :)
    eval_sdf!(o.Ω, c, points)
    c .+= o.lim

    # return obj - o.μ * sum(c -> log(-c), c)
    return obj + o.μ * sum(c -> 1 / (abs(c)), c)
end

function opt_gradient!(o::BFGSOptimizerState, g, x, c)
    # Clear gradient
    g .= 0.0

    # Compute penalty gradient
    points = reshape(x, 3, :)
    g_points = reshape(g, 3, :)

    eval_sdf!(o.Ω, c, points)
    # @show c
    c .+= o.lim
    eval_sdf_grad!(o.Ω, g_points, points)

    g_points .*= o.μ ./ c'
    g_points ./= c' # For 1/c

    # Compute main objective gradient
    o.g!(g, x)
end

function descent_direction!(o::BFGSOptimizerState; max_reg_iter=50, ρₘᵢₙ=1e-10, ρₘₐₓ=1e10)
    Bₖ = o.Bₖ
    Bₒ = o.Bₒ

    # Ensure B is symmetric/Hermitian
    Bₖ .+= Bₖ'
    Bₖ ./= 2.0

    # Back up Hessian estimate
    Bₒ .= Bₖ

    for k in 1:max_reg_iter
        # Add in regularization
        Bₖ[diagind(Bₖ)] .+= o.ρ

        # Perform Cholesky factorization, don't check for positive definiteness
        B_factorized = LinearAlgebra.cholesky!(Hermitian(Bₖ), check=false)

        if issuccess(B_factorized)
            o.ρ = clamp(0.8 * o.ρ, 0.0, ρₘₐₓ)

            # Use factorization to compute descent direction
            ldiv!(o.p, B_factorized, o.gₖ)
            o.p .*= -1.0

            break
        else
            Bₖ .= Bₒ
            o.ρ = clamp(1.6 * o.ρ, ρₘᵢₙ, ρₘₐₓ)
        end
    end

    # Restore pre-regularized Hessian estimate
    Bₖ .= Bₒ
end


function linesearch!(o::BFGSOptimizerState)
    o.α = clamp(1.5 * o.α, 0.0, 1.0)

    while o.α ≥ 1e-12
        @. o.xₖ₊₁ = o.xₖ + o.α * o.p
        points = reshape(o.xₖ₊₁, 3, :)
        eval_sdf!(o.Ω, o.cₖ₊₁, points)

        if (maximum(o.cₖ₊₁) < -o.lim)
            o.fₖ₊₁ = opt_objective(o, o.xₖ₊₁, o.cₖ₊₁)

            if (o.fₖ₊₁ ≤ o.fₖ + 1e-4 * o.α * dot(o.p, o.gₖ))
                @. o.s = o.α * o.p
                return true
                # else
                #     println("insufficient descent")
            end
            # else
            #     println("overshoot")
        end

        o.α *= 0.2
    end
    println("Line-search failed, α: $(o.α)")
    o.α = 1e-13
    @. o.xₖ₊₁ = o.xₖ + o.α * o.p
    @. o.s = o.α * o.p
    o.fₖ₊₁ = opt_objective(o, o.xₖ₊₁, o.cₖ₊₁)

    o.α = 1.0
    return false
end

# BFGS optimizer for voronoi tessellations, with line-search to ensure strict feasibility
function optimize_voronoi!(o::BFGSOptimizerState, x₀)
    # Initialize from initial guess
    initialize!(o, x₀)

    for k in 1:o.max_iters
        # Get gradient at xₖ from previous iteration
        o.xₖ .= o.xₖ₊₁
        o.gₖ .= o.gₖ₊₁
        o.cₖ .= o.cₖ₊₁
        o.fₖ = o.fₖ₊₁

        # Solve for descent direction pₖ
        descent_direction!(o)

        # Perform line-search to ensure feasibility and sufficient descent
        if !linesearch!(o)
            # println("Resetting hessian")
            reset_hessian!(o)
            opt_gradient!(o, o.gₖ, o.xₖ, o.cₖ)
            descent_direction!(o)
            linesearch!(o)
        end

        # Iterate updated by linesearch

        # Get gradient at xₖ₊₁, compute secant
        opt_gradient!(o, o.gₖ₊₁, o.xₖ₊₁, o.cₖ₊₁)
        @. o.y = o.gₖ₊₁ - o.gₖ

        # sTy = dot(o.s, o.y)
        sBs = dot(o.s' * o.Bₖ, o.s)

        # Compute curvature damping
        # θ = (sTy ≥ 0.2 * sBs) ? 1.0 : (0.8 * sBs) / (sBs - sTy)

        # Update damped secant
        # r = θ * y + (1 - θ) * B * s
        # mul!(o.r, o.Bₖ, o.s)
        # o.r .*= (1.0 - θ)
        # o.r .+= θ .* o.y
        o.r .= o.y

        sTr = dot(o.s, o.r)

        # Update Hessian estimate
        # B = B - (B * s * s' * B) / (s' * B * s) + (r * r') / (s' * r)
        o.Bₒ .= o.Bₖ
        mul!(o.o_prod_1, o.s, o.s')
        mul!(o.o_prod_2, o.Bₒ, o.o_prod_1)
        mul!(o.Bₖ, o.o_prod_2, o.Bₒ, -1.0 / sBs, 1.0)
        mul!(o.Bₖ, o.r, o.r', 1.0 / sTr, 1.0)

        # Print out iteration update
        # @printf "k    f         logmu α  ρ\n"
        @printf "%4i % 13.6e %.3e %.3e %.3e %.3e\n" k o.fₖ₊₁ o.μ o.α o.ρ norm(o.gₖ₊₁)

        if abs(o.fₖ - o.fₖ₊₁) < o.f_abs_tol
            println("Absolute function tolerance satisfied.")
            if o.μ > 1e-6
                o.μ *= 0.3
                o.fₖ₊₁ = opt_objective(o, o.xₖ₊₁, o.cₖ₊₁)
                opt_gradient!(o, o.gₖ₊₁, o.xₖ₊₁, o.cₖ₊₁)
                reset_hessian!(o)
                continue
            else
                println("Optimal solution found")
                break
            end
        end

        if abs((o.fₖ - o.fₖ₊₁) / o.fₖ) < o.f_rel_tol
            println("Relative function tolerance satisfied.")
            if o.μ > 1e-6
                o.μ *= 0.3
                o.fₖ₊₁ = opt_objective(o, o.xₖ₊₁, o.cₖ₊₁)
                opt_gradient!(o, o.gₖ₊₁, o.xₖ₊₁, o.cₖ₊₁)
                reset_hessian!(o)
                continue
            else
                println("Optimal solution found")
                break
            end
        end

        if norm(o.gₖ) / sqrt(length(o.gₖ)) < o.g_abs_tol
            println("Absolute gradient tolerance satisfied.")
            if o.μ > 1e-6
                o.μ *= 0.3
                o.fₖ₊₁ = opt_objective(o, o.xₖ₊₁, o.cₖ₊₁)
                opt_gradient!(o, o.gₖ₊₁, o.xₖ₊₁, o.cₖ₊₁)
                reset_hessian!(o)
                continue
            else
                println("Optimal solution found")
                break
            end
        end

        # @show o.gₖ
        # if o.gₖ' * inv(o.Bₖ) * o.gₖ < o.g_rel_tol
        #     println("Relative gradient tolerance satisfied.")
        #     if o.μ > 1e-6
        #         o.μ *= 0.3
        #         o.fₖ₊₁ = opt_objective(o, o.xₖ₊₁, o.cₖ₊₁)
        #         opt_gradient!(o, o.gₖ₊₁, o.xₖ₊₁, o.cₖ₊₁)
        #         reset_hessian!(o)
        #         continue
        #     else
        #         println("Optimal solution found")
        #         break
        #     end
        # end
    end
end
