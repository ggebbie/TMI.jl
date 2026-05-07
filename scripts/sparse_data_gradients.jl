import Pkg; Pkg.activate(".")

using TMI
using Enzyme
using Interpolations
using LinearAlgebra
using Printf
using BenchmarkTools

"""
Compare manual adjoint gradients and Enzyme reverse-mode gradients for
`costfunction_point_obs!` on sparse (pointwise) observations.
"""

"""
    point_obs_cost_with_grad!(..., guvec)

Scalar cost wrapper that writes `∂J/∂uvec` into `guvec`.
"""
function point_obs_cost_with_grad!(
    uvec::Vector, Alu, b, u₀, y::Vector, W⁻::Diagonal, wis, locs, Q⁻, γ::Grid, guvec,
)
    return TMI.costfunction_point_obs!(0.0, guvec, uvec, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ)
end

"""
    manual_gradient!(grad, ...)

Reference gradient path using TMI's hand-written adjoint.
"""
function manual_gradient!(
    grad, uvec::Vector, Alu, b, u₀, y::Vector, W⁻::Diagonal, wis, locs, Q⁻, γ::Grid,
)
    fill!(grad, 0)
    J = point_obs_cost_with_grad!(uvec, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ, grad)
    return J, grad
end

"""
    enzyme_gradient!(grad, ...)

Enzyme reverse-mode gradient path for the same scalar objective.
Only `uvec` is active; model/data objects are constants.
"""
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

manual_trial = @benchmark manual_wrapper!($uvec, $manual_grad)
enzyme_trial = @benchmark enzyme_wrapper!($uvec, $grad_uvec_ad)

manual_time = median(manual_trial).time * 1e-9
enzyme_time = median(enzyme_trial).time * 1e-9

absolute_error = maximum(abs.(grad_uvec_ad .- manual_grad))
relative_error = absolute_error / max(norm(manual_grad, Inf), eps())

println()
println("Enzyme vs. manual gradient comparison for costfunction_point_obs!")
println("---------------------------------------------------")
@printf("Objective J:                   %.6e\n", J)
@printf("Max absolute grad error:       %.6e\n", absolute_error)
@printf("Max relative grad error:       %.6e\n", relative_error)
println()
@printf("Manual gradient runtime:       %.3f ms\n", manual_time * 1e3)
@printf("Enzyme gradient runtime:       %.3f ms\n", enzyme_time * 1e3)
@printf("Runtime ratio (Enzyme/manual): %.2fx\n", enzyme_time / manual_time)
println("---------------------------------------------------")
