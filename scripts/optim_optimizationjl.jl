import Pkg; Pkg.activate(".")

using Optimization
using OptimizationOptimJL
using LinearAlgebra
using LineSearches
using Optim

# Extended Powell function
function powell(x)
    n = length(x)
    @assert n % 4 == 0
    s = 0.0
    for i in 1:4:n
        s += (x[i] + 10x[i+1])^2
        s += 5 * (x[i+2] - x[i+3])^2
        s += (x[i+1] - 2x[i+2])^4
        s += 10 * (x[i] - x[i+3])^4
    end
    return s
end

# Gradient of Powell function
function powell_grad!(g, x)
    n = length(x)
    @assert n % 4 == 0
    g .= 0.0
    for i in 1:4:n
        g[i]     += 2 * (x[i] + 10x[i+1]) + 40 * (x[i] - x[i+3])^3
        g[i+1]   += 20 * (x[i] + 10x[i+1]) + 4 * (x[i+1] - 2x[i+2])^3
        g[i+2]   += 10 * (x[i+2] - x[i+3]) - 8 * (x[i+1] - 2x[i+2])^3
        g[i+3]   += -10 * (x[i+2] - x[i+3]) - 40 * (x[i] - x[i+3])^3
    end
    return g
end

# Problem size and bounds
n = 4^5  # must be divisible by 4
x0 = abs.(rand(n))              # Ensure strictly inside bounds
lower = fill(0.0, n)
upper = fill(10.0, n)

# Build OptimizationFunction with separate f and g; Optimization.jl will use both.
optf = OptimizationFunction((x, p) -> powell(x); grad =  (G, x, p) -> powell_grad!(G, x))
prob = OptimizationProblem(optf, x0, nothing; lb = lower, ub = upper)

# Use Optim.jl's bounded LBFGS via the OptimJL backend (options can be added later)
sol = Optimization.solve(prob, OptimizationOptimJL.Fminbox(OptimizationOptimJL.LBFGS()))

# Alternative: single fg! following the Optim.jl template.
function powell_fg_optim!(F, G, x)
    fval = powell(x)
    if G !== nothing
        powell_grad!(G, x)
    end
    if F !== nothing
        return fval
    end
    return fval
end

result_opt_fg = Optim.optimize(
    Optim.only_fg!(powell_fg_optim!),
    lower,
    upper,
    x0,
    Optim.Fminbox(Optim.LBFGS())
)

# Compare to calling Optim.jl directly for reference.
od = Optim.OnceDifferentiable(powell, powell_grad!, x0)
result_opt = Optim.optimize(od, lower, upper, x0,
    Optim.Fminbox(Optim.LBFGS()), 
        Optim.Options(
        iterations = 2500,
        store_trace = false,
        show_trace = false,
        show_warnings = false))



# Direct Optim.jl with custom options mirroring scripts/optim.jl
result_opt_opts = Optim.optimize(
    od, lower, upper, x0,
    Optim.Fminbox(Optim.LBFGS(; m = 10,
        alphaguess = LineSearches.InitialHagerZhang(α0 = NaN),
        linesearch = LineSearches.HagerZhang())),
    Optim.Options(
        g_tol = 1e-12,
        outer_iterations = 1000,
        iterations = 2000,
        store_trace = false,
        show_trace = false,
        show_warnings = false)
)

println("Optimization.jl + OptimJL result:")
println("  f(x) = ", sol.objective)
println("  ||x|| = ", norm(sol.u))

println("\nOptim.jl result (fg!):")
println("  f(x) = ", result_opt_fg.minimum)
println("  ||x|| = ", norm(result_opt_fg.minimizer))

println("\nDirect Optim.jl result:")
println("  f(x) = ", result_opt.minimum)
println("  ||x|| = ", norm(result_opt.minimizer))

println("\nDirect Optim.jl with custom options:")
println("  f(x) = ", result_opt_opts.minimum)
println("  ||x|| = ", norm(result_opt_opts.minimizer))
