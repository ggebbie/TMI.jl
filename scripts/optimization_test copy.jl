# Benchmarking IPOPT vs Optim.Fminbox(LBFGS) on high-dimensional Powell function

using Optimization, OptimizationIpopt, OptimizationOptimJL
using LinearAlgebra, Statistics, Random, Printf
using BenchmarkTools, Plots

# Extended Powell function and gradient for n divisible by 4
function powell(x, p)
    n = length(x)
    @assert n % 4 == 0 "Powell function requires n divisible by 4"
    s = 0.0
    for i in 1:4:n
        s += (x[i] + 10x[i+1])^2
        s += 5 * (x[i+2] - x[i+3])^2
        s += (x[i+1] - 2x[i+2])^4
        s += 10 * (x[i] - x[i+3])^4
    end
    return s
end

function powell_grad!(G, x, p)
    n = length(x)
    @assert n % 4 == 0
    G .= 0.0
    for i in 1:4:n
        # term1: (x[i] + 10x[i+1])^2
        G[i]     += 2 * (x[i] + 10x[i+1])
        G[i+1]   += 20 * (x[i] + 10x[i+1])

        # term2: 5(x[i+2] - x[i+3])^2
        G[i+2]   += 10 * (x[i+2] - x[i+3])
        G[i+3]   += -10 * (x[i+2] - x[i+3])

        # term3: (x[i+1] - 2x[i+2])^4
        t3 = (x[i+1] - 2x[i+2])^3
        G[i+1]   += 4 * t3
        G[i+2]   += -8 * t3

        # term4: 10(x[i] - x[i+3])^4
        t4 = (x[i] - x[i+3])^3
        G[i]     += 40 * t4
        G[i+3]   += -40 * t4
    end
end

# Benchmark over increasing problem sizes
ns = [4^4, 4^5, 4^6, 4^7]  # Must be divisible by 4
ipopt_times = zeros(length(ns))
optim_times = zeros(length(ns))
ipopt_errors = zeros(length(ns))
optim_errors = zeros(length(ns))
ipopt_memory = zeros(length(ns))
optim_memory = zeros(length(ns))

x0 = abs.(randn(n))
n = 4^5
lb = fill(0.0, n)
ub = fill(10.0, n)
p = nothing  # unused, for API compatibility
true_sol = zeros(n)  # Known minimizer of Powell function

optfun = OptimizationFunction(powell; grad = powell_grad!)
prob = OptimizationProblem(optfun, x0, p; lb = lb, ub = ub)
optim_bench = solve(prob, Optim.Fminbox(Optim.LBFGS(m=10); Optim.Options(outer_iterations = 2)); maxiters = 1000)
    
println("\nBenchmarking IPOPT vs Optim.Fminbox(LBFGS) on Extended Powell Function")

for (j, n) in enumerate(ns)
    println("\n--- n = $n ---")
    x0 = randn(n)
    lb = fill(-10.0, n)
    ub = fill(10.0, n)
    p = nothing  # unused, for API compatibility
    true_sol = zeros(n)  # Known minimizer of Powell function

    optfun = OptimizationFunction(powell; grad = powell_grad!)
    prob = OptimizationProblem(optfun, x0, p; lb = lb, ub = ub)

    # IPOPT benchmark
    println("Running IPOPT...")
    ipopt_bench = @benchmark begin
        solve($prob, IpoptOptimizer(
            hessian_approximation = "limited-memory",
            limited_memory_max_history = 10);
            max_iter = 1000,
            tol = 1e-6,
            print_level = 0
        )
    end samples=1 evals=1
    ipopt_times[j] = ipopt_bench.times[1] / 1e9
    ipopt_memory[j] = ipopt_bench.memory / 1e9  # Convert to GB
    ipopt_sol = solve(prob, IpoptOptimizer(
        hessian_approximation = "limited-memory",
        limited_memory_max_history = 10);
        max_iter = 1000,
        tol = 1e-6,
        print_level = 0
    )
    ipopt_errors[j] = norm(ipopt_sol.u .- true_sol)

    optfun = OptimizationFunction(powell; grad = powell_grad!)
    prob = OptimizationProblem(optfun, x0, p; lb = lb, ub = ub)

    # Optim benchmark
    println("Running Optim.Fminbox(LBFGS)...")
    optim_bench = @benchmark begin
        solve($prob, Optim.Fminbox(Optim.LBFGS()); maxiters = 1000)
    end samples=1 evals=1
    optim_times[j] = optim_bench.times[1] / 1e9
    optim_memory[j] = optim_bench.memory / 1e9  # Convert to GB
    optim_sol = solve(prob, Optim.Fminbox(Optim.LBFGS()); maxiters = 1000)
    optim_errors[j] = norm(optim_sol.u .- true_sol)
end

# Subplots for timing, error, and memory
# plt1 = plot(ns, ipopt_times, label = "IPOPT", lw = 2, marker = :circle);
plt1 = plot(ns, optim_times, label = "Fminbox(LBFGS)", lw = 2, marker = :diamond);
xlabel!(plt1, "Problem Dimension (n)");
ylabel!(plt1, "Solve Time (s)");
title!(plt1, "Extended Powell Benchmark: Timing");

# plt2 = plot(ns, ipopt_errors, label = "IPOPT Error", lw = 2, marker = :circle);
plt2 = plot(ns, optim_errors, label = "Fminbox Error", lw = 2, marker = :diamond);
xlabel!(plt2, "Problem Dimension (n)");
ylabel!(plt2, "||sol - sol_true||");
title!(plt2, "Distance to True Solution");

# plt3 = plot(ns, ipopt_memory, label = "IPOPT", lw = 2, marker = :circle);
plt3 = plot(ns, optim_memory, label = "Fminbox(LBFGS)", lw = 2, marker = :diamond);
xlabel!(plt3, "Problem Dimension (n)");
ylabel!(plt3, "Memory (GB)");
title!(plt3, "Memory Usage");
# grid!(plt3, true);

plot(plt1, plt2, plt3, layout = (1, 3), size = (1400, 400))


# Subplots for timing, error, and memory
plt1 = plot(ns, ipopt_times, label = "IPOPT", lw = 2, marker = :circle);
plot!(plt1, ns, optim_times, label = "Fminbox(LBFGS)", lw = 2, marker = :diamond);
xlabel!(plt1, "Problem Dimension (n)");
ylabel!(plt1, "Solve Time (s)");
title!(plt1, "Extended Powell Benchmark: Timing");

plt2 = plot(ns, ipopt_errors, label = "IPOPT Error", lw = 2, marker = :circle);
plot!(plt2, ns, optim_errors, label = "Fminbox Error", lw = 2, marker = :diamond);
xlabel!(plt2, "Problem Dimension (n)");
ylabel!(plt2, "||sol - sol_true||");
title!(plt2, "Distance to True Solution");

plt3 = plot(ns, ipopt_memory, label = "IPOPT", lw = 2, marker = :circle);
plot!(plt3, ns, optim_memory, label = "Fminbox(LBFGS)", lw = 2, marker = :diamond);
xlabel!(plt3, "Problem Dimension (n)");
ylabel!(plt3, "Memory (GB)");
title!(plt3, "Memory Usage");
# grid!(plt3, true);

plot(plt1, plt2, plt3, layout = (1, 3), size = (1400, 400))
