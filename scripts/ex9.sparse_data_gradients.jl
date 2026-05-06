import Pkg; Pkg.activate(".")

using TMI
using Enzyme
using LinearAlgebra
using BenchmarkTools

function costfunction_point_obs_value(uvec::Vector, Alu, b, u₀, y::Vector, W⁻::Diagonal, locs, Q⁻, γ::Grid, guvec)
    return TMI.costfunction_point_obs!(0.0, guvec, uvec, Alu, b, u₀, y, W⁻, locs, locs, Q⁻, γ)
end

function enzyme_gradient!(grad, p)
    fill!(grad, 0)
    Enzyme.autodiff(set_runtime_activity(Reverse), costfunction_point_obs_value, Duplicated(p.uvec, grad), Const(p.Alu), Const(p.b), Const(p.u₀), Const(p.y), Const(p.W⁻), Const(p.locs), Const(p.Q⁻), Const(p.γ), Const(nothing))
    return grad
end

function manual_gradient(p)
    guvec = zero(p.uvec)
    J = costfunction_point_obs_value(p.uvec, p.Alu, p.b, p.u₀, p.y, p.W⁻, p.locs, p.Q⁻, p.γ, guvec)
    return J, guvec
end

TMIversion = "modern_90x45x33_GH10_GH12"
N = 20
σb = 5.0
correlation_length = 1000.0

A, Alu, γ, TMIfile, L, B = config(TMIversion)

u₀ = (;surface = zerosurfaceboundary(γ))
# u₀  = zerosurfaceboundary(γ)

uvec = copy(vec(u₀))

y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion, "θ", γ, N)
b = (;surface = mean(y) * onesurfaceboundary(γ))
# b  = mean(y) * onesurfaceboundary(γ)

Dg = gaussiandistancematrix(γ, σb, correlation_length)
Q⁻ = inv(cholesky(Dg))

parameters = (; uvec, Alu, b, u₀, y, W⁻, wis, locs, Q⁻, γ)

J, guvec = manual_gradient(parameters)

grad_uvec_ad = Enzyme.make_zero(uvec)
enzyme_gradient!(grad_uvec_ad, parameters)

absolute_error = maximum(abs.(grad_uvec_ad .- guvec))
manual_benchmark = @benchmark manual_gradient($parameters)
enzyme_benchmark = @benchmark enzyme_gradient!($grad_uvec_ad, $parameters)

println("Gradient benchmark summary:")
println("  J = ", J)
println("  absolute_error = ", absolute_error)
println()
println("Manual gradient:")
display(manual_benchmark)
println()
println("Enzyme gradient:")
display(enzyme_benchmark)
