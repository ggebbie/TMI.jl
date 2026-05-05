import Pkg; Pkg.activate(".")

using TMI
using BenchmarkTools
using Enzyme
using LinearAlgebra
using Statistics

function add_boundarycondition_value(b::BoundaryCondition, u::BoundaryCondition)
    tracer = b.tracer .+ u.tracer
    return BoundaryCondition(tracer, b.axes, b.k, b.dim, b.dimval, b.wet, b.name, b.longname, b.units)
end

function costfunction_point_obs_value(uvec::Vector, Alu, b::BoundaryCondition, u₀::BoundaryCondition, y::Vector, Wⁱ::Diagonal, locs, Q⁻, γ::Grid)
    u = unvec(u₀, uvec)
    b̃ = add_boundarycondition_value(b, u)
    c = steadyinversion(Alu, b̃, γ)
    ỹ = observe(c, locs, γ)
    n = ỹ - y
    return n ⋅ (Wⁱ * n) + uvec' * (Q⁻ * uvec)
end

function enzyme_gradient!(grad, p)
    fill!(grad, 0)
    args = (
        Duplicated(p.uvec, grad), Const(p.Alu), Const(p.b), Const(p.u),
        Const(p.y), Const(p.W⁻), Const(p.locs), Const(p.Q⁻), Const(p.γ),
    )
    Enzyme.autodiff(Enzyme.set_runtime_activity(Reverse), costfunction_point_obs_value, args...)
    return grad
end

manual_gradient(p) = TMI.costfunction_point_obs(p.uvec, p.Alu, p.b, p.u, p.y, p.W⁻, p.wis, p.locs, p.Q⁻, p.γ)

TMIversion = "modern_90x45x33_GH10_GH12"
N = 20
σb = 5.0
correlation_length = 1000.0

A, Alu, γ, TMIfile, L, B = config(TMIversion)

u = zerosurfaceboundary(γ)
uvec = vec(u)

y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion, "θ", γ, N)
b = mean(y) * onesurfaceboundary(γ)

Dg = gaussiandistancematrix(γ, σb, correlation_length)
Q⁻ = inv(cholesky(Dg))

problem = (; uvec, Alu, b, u, y, W⁻, wis, locs, Q⁻, γ)

J, guvec = manual_gradient(problem)

grad_b_ad = Enzyme.make_zero(uvec)
enzyme_gradient!(grad_b_ad, problem)
absolute_error = maximum(abs.(grad_b_ad .- guvec))

manual_time = @belapsed manual_gradient($problem) samples=3 evals=1
enzyme_time = @belapsed enzyme_gradient!($grad_b_ad, $problem) samples=3 evals=1

println("Manual gradient benchmark:")
@btime manual_gradient($problem) samples=3 evals=1

println("Enzyme gradient benchmark:")
@btime enzyme_gradient!($grad_b_ad, $problem) samples=3 evals=1

println("Summary:")
println("  J = ", J)
println("  manual_time = ", manual_time)
println("  enzyme_time = ", enzyme_time)
println("  absolute_error = ", absolute_error)
