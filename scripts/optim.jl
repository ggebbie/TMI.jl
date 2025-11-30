using Optim, BenchmarkTools, LinearAlgebra
using LineSearches
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

# Problem size
n = 4^5  # must be divisible by 4
x0 = abs.(rand(n))              # Ensure strictly inside bounds
lower = fill(0.0, n)
upper = fill(10.0, n)

# Set up OnceDifferentiable
od = OnceDifferentiable(powell, powell_grad!, x0)

# Run optimization
# @btime result = Optim.optimize($od, $lower, $upper, $x0, Fminbox(LBFGS()));
result1 = Optim.optimize(od, lower, upper, x0, 
                            Fminbox(LBFGS(;m = 10, 
                              alphaguess = LineSearches.InitialHagerZhang(α0=NaN),
                                linesearch = LineSearches.HagerZhang())),
                            Optim.Options(
                            f_tol = 1e-12,
                            g_tol = 1e-12,
                             iterations = 2500,
                             store_trace = false,
                             show_trace = false,
                             show_warnings = false));
# Outpprintln("\nFinal Objective Value: ", result.minimum)
println("||sol - sol_true|| = ", norm(result1.minimizer))  # True solution is 0 vector
result1