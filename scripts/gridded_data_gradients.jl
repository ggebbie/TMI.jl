import Pkg; Pkg.activate(".")

using TMI
using Enzyme
using LinearAlgebra
using Printf
using BenchmarkTools

"""
Compare manual adjoint gradients and Enzyme reverse-mode gradients for
`costfunction_gridded_obs!` on gridded observations.
"""

"""
    gridded_obs_cost_with_grad!(..., c0_work, ..., guvec)

Scalar cost wrapper that writes `∂J/∂uvec` into `guvec`.
`c0_work` is mutated by `costfunction_gridded_obs!`.
"""
function gridded_obs_cost_with_grad!(
    uvec::Vector, Alu, b, u₀, c0_work::Field, W⁻::Diagonal, γ::Grid, guvec,
)
    return TMI.costfunction_gridded_obs!(0.0, guvec, uvec, Alu, b, u₀, c0_work, W⁻, γ)
end

"""
    manual_gradient!(grad, ...)

Reference gradient path using TMI's hand-written adjoint.
"""
function manual_gradient!(
    grad, uvec::Vector, Alu, b, u₀, c0::Field, W⁻::Diagonal, γ::Grid,
)
    fill!(grad, 0)
    c0_work = deepcopy(c0)
    J = gridded_obs_cost_with_grad!(uvec, Alu, b, u₀, c0_work, W⁻, γ, grad)
    return J, grad
end

"""
    enzyme_gradient!(grad, ...)

Enzyme reverse-mode gradient path for the same scalar objective.
Only `uvec` is active; model/data objects are constants.
"""
function enzyme_gradient!(
    grad, uvec::Vector, Alu, b, u₀, c0::Field, W⁻::Diagonal, γ::Grid,
)
    fill!(grad, 0)
    c0_work = deepcopy(c0)
    Enzyme.autodiff(
        set_runtime_activity(Reverse),
        gridded_obs_cost_with_grad!,
        Duplicated(uvec, grad),
        Const(Alu),
        Const(b),
        Const(u₀),
        Const(c0_work),
        Const(W⁻),
        Const(γ),
        Const(nothing),
    )
    return grad
end

TMIversion = "modern_90x45x33_GH10_GH12"

A, Alu, γ, TMIfile, L, B = config(TMIversion)

u₀ = (;surface = zerosurfaceboundary(γ))
uvec = copy(vec(u₀))

θobs = readfield(TMIfile, "θ", γ)
b_truth = (;surface = getsurfaceboundary(θobs))
c0 = steadyinversion(Alu, b_truth, γ)

σθ = readfield(TMIfile, "σθ", γ)
W⁻ = (1 / sum(γ.wet)) .* Diagonal(1 ./σθ.tracer[γ.wet].^2)

b = (;surface = zerosurfaceboundary(γ))

manual_grad = zero(uvec)
J, _ = manual_gradient!(manual_grad, uvec, Alu, b, u₀, c0, W⁻, γ)

grad_uvec_ad = Enzyme.make_zero(uvec)
enzyme_gradient!(grad_uvec_ad, uvec, Alu, b, u₀, c0, W⁻, γ)

manual_wrapper!(uvec_in, grad_out) = manual_gradient!(grad_out, uvec_in, Alu, b, u₀, c0, W⁻, γ)
enzyme_wrapper!(uvec_in, grad_out) = enzyme_gradient!(grad_out, uvec_in, Alu, b, u₀, c0, W⁻, γ)

manual_trial = @benchmark manual_wrapper!($uvec, $manual_grad)
enzyme_trial = @benchmark enzyme_wrapper!($uvec, $grad_uvec_ad)

manual_time = median(manual_trial).time * 1e-9
enzyme_time = median(enzyme_trial).time * 1e-9

absolute_error = maximum(abs.(grad_uvec_ad .- manual_grad))
relative_error = absolute_error / max(norm(manual_grad, Inf), eps())

println()
println("Enzyme vs. manual gradient comparison for costfunction_gridded_obs!")
println("---------------------------------------------------")
@printf("Objective J:              %.6e\n", J)
@printf("Max absolute grad error:  %.6e\n", absolute_error)
@printf("Max relative grad error:  %.6e\n", relative_error)
println()
@printf("Manual gradient:          %.3f ms\n", manual_time * 1e3)
@printf("Enzyme gradient:          %.3f ms\n", enzyme_time * 1e3)
@printf("Runtime ratio (Enzyme/manual): %.2fx\n", enzyme_time / manual_time)
println("---------------------------------------------------")
