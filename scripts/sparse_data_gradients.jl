import Pkg; Pkg.activate(".")

using TMI
using Enzyme
using Interpolations
using LinearAlgebra
using BenchmarkTools

function point_obs_cost_with_grad!(
    uvec::Vector, Alu, b, u₀, y::Vector, W⁻::Diagonal, wis, locs, Q⁻, γ::Grid, guvec,
)
    return TMI.costfunction_point_obs!(0.0, guvec, uvec, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ)
end

function manual_gradient!(
    grad, uvec::Vector, Alu, b, u₀, y::Vector, W⁻::Diagonal, wis, locs, Q⁻, γ::Grid,
)
    fill!(grad, 0)
    J = point_obs_cost_with_grad!(uvec, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ, grad)
    return J, grad
end

function enzyme_gradient!(
    grad, uvec::Vector, Alu, b, u₀, y::Vector, W⁻::Diagonal, wis, locs, Q⁻, γ::Grid,
)
    fill!(grad, 0)
    Enzyme.autodiff(
        set_runtime_activity(Reverse),
        point_obs_cost_with_grad!,
        Duplicated(uvec, grad),
        Const(Alu),
        Const(b),
        Const(u₀),
        Const(y),
        Const(W⁻),
        Const(wis),
        Const(locs),
        Const(Q⁻), Const(γ), Const(nothing),
    )
    return grad
end

TMIversion = "modern_90x45x33_GH10_GH12"
N = 20
σb = 5.0
correlation_length = 1000.0

A, Alu, γ, TMIfile, L, B = config(TMIversion)

u₀ = (;surface = zerosurfaceboundary(γ))

uvec = copy(vec(u₀))

y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion, "θ", γ, N)
b = (;surface = mean(y) * onesurfaceboundary(γ))

Dg = gaussiandistancematrix(γ, σb, correlation_length)
Q⁻ = inv(cholesky(Dg))

manual_grad = zero(uvec)
J, _ = manual_gradient!(manual_grad, uvec, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ)

grad_uvec_ad = Enzyme.make_zero(uvec)
enzyme_gradient!(grad_uvec_ad, uvec, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ)

absolute_error = maximum(abs.(grad_uvec_ad .- manual_grad))

manual_wrapper!(uvec_in, grad_out) = manual_gradient!(grad_out, uvec_in, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ)
enzyme_wrapper!(uvec_in, grad_out) = enzyme_gradient!(grad_out, uvec_in, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ)

println("Gradient benchmark summary:")
println("  J = ", J)
println("  absolute_error = ", absolute_error)
println()
println("Manual gradient (@btime):")
@btime manual_wrapper!($uvec, $manual_grad)
println()
println("Enzyme gradient (@btime):")
@btime enzyme_wrapper!($uvec, $grad_uvec_ad)
