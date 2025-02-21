function lloyd(points::DMT, Ω::F; SMT=SparseMatrixCSC{Int,Int}, tol=1e-3, max_iters=1000) where {DMT,F}
    n = size(points, 2)
    problem = CVT_Problem(Ω, n)

    solver = Ipopt.Optimizer()
    solver.options["max_iter"] = max_iters
    solver.options["tol"] = tol
    solver.options["constr_viol_tol"] = tol
    # solver.options["hsllib"] = HSL_jll.libhsl_path
    # solver.options["linear_solver"] = solv
    solver.options["neg_curv_test_tol"] = 1e-11

    # if isqp
    #     solver.options["hessian_constant"] = "yes"
    #     solver.options["jac_c_constant"] = "yes"
    #     solver.options["jac_d_constant"] = "yes"
    # end

    ct = zeros(n)
    gt = zeros(3n)
    jact = zeros(3n)
    z₀ = vec(points)

    # if isnothing(problem.hessian_lagrangian_sparsity)
    #     println("Building sparsity for the lagrangian hessian...")
    #     build_hessian_lagrangian_structure!(problem)
    # end

    println("Checking objective function...")
    @time MathOptInterface.eval_objective(problem, z₀)
    @time MathOptInterface.eval_objective(problem, z₀)

    println("Checking objective gradient...")
    @time MathOptInterface.eval_objective_gradient(problem, gt, z₀)
    @time MathOptInterface.eval_objective_gradient(problem, gt, z₀)
    @show gt

    println("Checking constraint function...")
    @time MathOptInterface.eval_constraint(problem, ct, z₀)
    @time MathOptInterface.eval_constraint(problem, ct, z₀)
    @show maximum(ct)

    println("Checking constraint jacobian...")
    @time MathOptInterface.eval_constraint_jacobian(problem, jact, z₀)
    @time MathOptInterface.eval_constraint_jacobian(problem, jact, z₀)
    @show maximum(jact)

    c_l = -Inf * ones(n)
    c_u = -0.1 * ones(n)

    nlp_bounds = MathOptInterface.NLPBoundsPair.(c_l, c_u)
    block_data = MathOptInterface.NLPBlockData(nlp_bounds, problem, true)

    z = MathOptInterface.add_variables(solver, 3n)

    # Set primal bounds and initial values
    # MathOptInterface.add_constraints(solver, z, MathOptInterface.LessThan.(100 * ones(3 * n)))
    # MathOptInterface.add_constraints(solver, z, MathOptInterface.GreaterThan.(zeros(3 * n)))
    # x_u[c][k] = MathOptInterface.add_constraints(solver, xj, MathOptInterface.LessThan.(problem.x̄))
    # x_l[c][k] = MathOptInterface.add_constraints(solver, xj, MathOptInterface.GreaterThan.(problem.xfmin))

    for i in 1:lastindex(z₀)
        MathOptInterface.set(solver, MathOptInterface.VariablePrimalStart(), z[i], z₀[i])
    end

    MathOptInterface.set(solver, MathOptInterface.NLPBlock(), block_data)
    MathOptInterface.set(solver, MathOptInterface.ObjectiveSense(), MathOptInterface.MIN_SENSE)

    flush(stdout)
    # Run optimizer
    MathOptInterface.optimize!(solver)

    # J = MathOptInterface.get(solver, MathOptInterface.ObjectiveValue())
    optimal_points = reshape(MathOptInterface.get(solver, MathOptInterface.VariablePrimal(), z), 3, :)
    # λ = MathOptInterface.get(solver, MathOptInterface.NLPBlockDual())
    # μ_uₗ[c][k] = MathOptInterface.get(solver, MathOptInterface.ConstraintDual(), u_l[c][k])

    return bounded_voronoi(optimal_points, Ω)
end

struct CVT_Problem{F} <: MathOptInterface.AbstractNLPEvaluator
    Ω::F
    n_points::Int
end

MathOptInterface.features_available(prob::CVT_Problem) = [:Grad, :Jac]#, :Hess]
MathOptInterface.initialize(prob::CVT_Problem, features) = nothing

function MathOptInterface.eval_objective(prob::CVT_Problem, z)
    points = reshape(z, 3, :)
    voronoi = bounded_voronoi(points, prob.Ω)

    J = 0.0
    for i in n_cells(voronoi)
        z_i = SVector{3}(points[1, i], points[2, i], points[3, i])
        f(y) = dot(y - z_i, y - z_i)
        J += cell_volume_integral(voronoi, f, i)
    end

    return J
end

function MathOptInterface.eval_objective_gradient(prob::CVT_Problem{F}, grad_f, z) where {F}
    points = reshape(z, 3, :)
    grad_p = reshape(grad_f, 3, :)

    voronoi = bounded_voronoi(points, prob.Ω)
    complex_centroids!(grad_p, voronoi)
    grad_p .*= -1
    grad_p .+= points
    @show any(isnan, grad_p)
    for i in 1:n_cells(voronoi)
        grad_p[:, i] .*= 2 * cell_volume(voronoi, i)
    end
end

function MathOptInterface.eval_constraint(prob::CVT_Problem, c, z)
    points = reshape(z, 3, :)
    map!(prob.Ω, c, eachslice(points, dims=2))
    c .*= 0.1
    @show maximum(c)
end

# Jacobian structure: row for each point, 3 columns for each point
function MathOptInterface.jacobian_structure(prob::CVT_Problem)
    cols = 1:(3*prob.n_points)
    rows = repeat(1:prob.n_points, inner=3)
    return collect(zip(rows, cols))
end

function MathOptInterface.eval_constraint_jacobian(prob::CVT_Problem{F}, jac, z) where {F}
    points = reshape(z, 3, :)
    jac_m = reshape(jac, 3, :)

    for i in 1:prob.n_points
        jac_i = view(jac_m, :, i)
        point_i = view(points, :, i)
        ForwardDiff.gradient!(jac_i, prob.Ω, point_i)
    end
end

# # Lagrangian Hessian structure: diagonal for the cost function, block diagonal for SDF
# function MathOptInterface.hessian_lagrangian_structure(prob::CVT_Problem)
#     rows_∇²f = 1:(3*prob.n_points)
#     cols_∇²f = 1:(3*prob.n_points)

#     rows_∂²c∂²z = vcat([repeat((3*(i-1)+1):(3*i), 3) for i in 1:prob.n_points]...)
#     cols_∂²c∂²z = vcat([repeat((3*(i-1)+1):(3*i), inner=3) for i in 1:prob.n_points]...)

#     rows = vcat(rows_∇²f, rows_∂²c∂²z)
#     cols = vcat(cols_∇²f, cols_∂²c∂²z)

#     return collect(zip(rows, cols))
# end

# function MathOptInterface.eval_hessian_lagrangian(prob::CVT_Problem{F}, H, z, σ, μ) where {F}
#     points = reshape(z, 3, :)
#     voronoi = bounded_voronoi(points, prob.Ω)

#     ∇²f = view(H, 1:(3*prob.n_points))
#     ∂²c∂²z = view(H, (3*prob.n_points+1):(length(H)))

#     ∇²f_m = reshape(∇²f, 3, :)
#     for i in 1:prob.n_points
#         ∇²f_m[:, i] .= 2 * cell_volume(voronoi, i)
#     end

#     ∇²f .*= σ

#     for i in 1:prob.n_points
#         h_i = view(∂²c∂²z, (9*(i-1)+1):9*i)
#         H_i = reshape(h_i, 3, 3)
#         point_i = view(points, :, i)

#         ForwardDiff.hessian!(H_i, prob.Ω, point_i)
#         H_i .*= μ[i]
#     end
# end
