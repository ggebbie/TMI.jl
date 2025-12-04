#!/usr/bin/env julia

using LinearAlgebra
using Optim
using NLopt
using BenchmarkTools

println("\n======== LARGE-SCALE (N=1000) FUSED-GRADIENT TEST ========\n")

############################################################
# Large-scale nonconvex problem (N = 1000)
#
#   f(x) = (x[i] - sin(i/50))^2
#        + 3*exp(-(x[i] - 0.8*sin(i/30))^2 / σ)
#
# Everything computed INSIDE fg_large! (no captured globals)
############################################################

const N = 1000
σ = 0.05

lb = fill(-5.0, N)
ub = fill( 5.0, N)

x0 = rand(N) .* 4 .- 2       # random initial guess in [-2,2]
@assert all(lb .<= x0 .<= ub)

println("Dimension: N = $N\n")


############################################################
#  FUSED GRADIENT CALLBACK (Anthony-style syntax)
############################################################

function fg_large!(F, G, x)
    # GRADIENT
    if G !== nothing
        # quadratic part
        for i in 1:length(x)
            target_i = sin(i/50)
            G[i] = 2*(x[i] - target_i)
        end

        # Gaussian wells
        for i in 1:length(x)
            well_i = 0.8*sin(i/30)
            dx = x[i] - well_i
            bump = exp(-(dx*dx)/σ)
            G[i] += bump * (-2/σ) * dx
        end
    end

    # OBJECTIVE
    if F !== nothing
        acc = 0.0
        for i in 1:length(x)
            target_i = sin(i/50)
            well_i   = 0.8*sin(i/30)
            dx1 = x[i] - target_i
            dx2 = x[i] - well_i
            acc += dx1*dx1 + 3*exp(-(dx2*dx2)/σ)
        end
        F = acc
        return acc
    end

    return nothing
end


############################################################
#  Optim.jl solver (LBFGS + Fminbox)
############################################################

function solve_optim(x)
    result = Optim.optimize(
        Optim.only_fg!(fg_large!),
        lb,
        ub,
        x,
        Fminbox(LBFGS()),
        Optim.Options(show_trace=false)
    )
    xmin = Optim.minimizer(result)
    fmin = Optim.minimum(result)
    return xmin, fmin
end

println("-- Optim.jl LBFGS --")
xmin_opt, fmin_opt = solve_optim(x0)
println("Optim minimizer:           ", xmin_opt)
println("Optim minimum objective:   ", fmin_opt, "\n")

println("Benchmarking Optim.jl...")
@btime solve_optim($x0);


############################################################
#  NLopt solver (LD_LBFGS)
############################################################

function solve_nlopt(x)
    opt = Opt(:LD_LBFGS, N)
    lower_bounds!(opt, lb)
    upper_bounds!(opt, ub)

    function objgrad(x, grad)
        # gradient
        if !isempty(grad)
            # quadratic
            for i in 1:N
                grad[i] = 2*(x[i] - sin(i/50))
            end
            # wells
            for i in 1:N
                well_i = 0.8*sin(i/30)
                dx = x[i] - well_i
                bump = exp(-(dx*dx)/σ)
                grad[i] += bump * (-2/σ) * dx
            end
        end

        # objective
        acc = 0.0
        for i in 1:N
            target_i = sin(i/50)
            well_i   = 0.8*sin(i/30)
            dx1 = x[i] - target_i
            dx2 = x[i] - well_i
            acc += dx1*dx1 + 3*exp(-(dx2*dx2)/σ)
        end
        return acc
    end

    min_objective!(opt, objgrad)
    xt = copy(x)
    (fmin, xmin, ret) = NLopt.optimize(opt, xt)
    return xmin, fmin
end

println("\n-- NLopt LD_LBFGS --")
xmin_nl, fmin_nl = solve_nlopt(x0)
println("NLopt minimizer:           ", xmin_nl)
println("NLopt minimum objective:   ", fmin_nl, "\n")

println("Benchmarking NLopt...")
@btime solve_nlopt($x0);

println("\n======== LARGE-SCALE TEST COMPLETE ========\n")
