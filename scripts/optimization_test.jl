import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using FiniteDiff
using Random
using Test

ngrid = (50) # number of grid cells
xmax = 1000.0 # domain size 
lon = collect(range(0.0,1000.0,length=ngrid[1]))
tracer = collect(1.0.-lon./xmax)

axes = (lon,)
wet = trues(ngrid)
interior = copy(wet)
interior[begin] = false
interior[end] = false

wrap = (false,)
Δ = [CartesianIndex(1,),CartesianIndex(-1,)]
γ = Grid(axes,wet,interior,wrap,Δ)
n = neighbors(γ)
m0 = massfractions_isotropic(γ)
m0 = (west = m0[1], east = m0[2])
c = Field(tracer,
    γ,
    :c,
    "linear equilibrated tracer",
    "μmol/kg")

A = watermassmatrix(m0, γ)
Alu = lu(A)

dim = 1
b = (west = TMI.getboundarycondition(c, dim, 1, γ),
        east = TMI.getboundarycondition(c, 1, ngrid[dim], γ))

# b = west = TMI.getboundarycondition(c, dim, 1, γ)
c̃ = steadyinversion(A,b,γ)
# Introduce a nonconservative variable with a source
qfield = 3.0e-2 * ones(ngrid)
# requires negative sign which is counterintutive (needs to be fixed)
q = TMI.Source(-qfield, γ, :q, "remineralized stuff", "μmol/kg", false)
c_noncons = steadyinversion(A,b,γ; q = q)
Δc = c_noncons - c

#deepcopy is necessary to avoid overalpping here 
u₀ = (c = deepcopy(b), c_q = deepcopy(b))
q₀ = (c = nothing, c_q = deepcopy(q),)
c0 = steadyinversion(Alu,u₀, q₀, γ)
Alu = lu(A)

# Build a single-point observation for each tracer using `observe`.
obs_loc = [lon[cld(length(lon), 2)]]
observed_vals = (; (k => observe(v, obs_loc, γ) for (k,v) in pairs(c0))...)
W = (; (k => Diagonal(fill!(similar(observed_vals[k], Float64), 1.0)) for k in keys(observed_vals))...)
c_obs = (; (k => Observations(observed_vals[k]; locs = obs_loc, γ = γ, W = W[k]) for k in keys(observed_vals))...)

du = zero(u₀)
dq = zero(q₀)

diag_from_v(v) = isnothing(v) ? nothing : Diagonal(one.(vec(v)))
Qᵤ = (; (k => diag_from_v(v) for (k,v) in pairs(du))...)
Qₛ = (; (k => diag_from_v(v) for (k,v) in pairs(dq))...)
Qₘ = Diagonal(one.(vec(m0)))

# Pull softmax helpers directly into this script for testing.
include(TMI.pkgsrcdir("parameterizations.jl"))

# Compare analytic adjoint of softmax to FiniteDiff on a simple squared cost.
function check_softmax_gradient(; n = 6, seed = 1)
    Random.seed!(seed)
    v = randn(n)

    softmax_cost(x) = sum(softmax_vector(x).^2)
    fd_grad = FiniteDiff.finite_difference_gradient(softmax_cost, v, Val{:central})

    s = softmax_vector(v)
    gs = 2 .* s # d/ds of sum(s.^2)
    analytic_grad = gsoftmax_vector(gs, s)

    println("softmax cost: ", softmax_cost(v))
    println("FiniteDiff gradient: ", fd_grad)
    println("analytic gradient: ", analytic_grad)
    @test fd_grad ≈ analytic_grad atol = 1e-10 rtol = 1e-10
    return nothing
end

check_softmax_gradient()

# MassFraction softmax: compare analytic adjoint to FiniteDiff using m0.
function check_massfraction_softmax_gradient(; seed = 3)
    Random.seed!(seed)
    mvec = randn(length(vec(m0)))

    cost_fn(x) = begin
        m = unvec(m0, x)
        sm = softmax_massfractions(m)
        # cost: sum of squared softmaxed mass fractions on interior points only
        return sum(sum(abs2, mf.fraction[mf.γ.interior]) for mf in sm)
    end

    fd_grad = FiniteDiff.finite_difference_gradient(cost_fn, mvec, Val{:central})

    # analytic gradient via adjoint
    m = unvec(m0, mvec)
    sm = softmax_massfractions(m)
    gs = map(sm) do mf
        gmf = similar(mf)
        gmf.fraction .= 0
        gmf.fraction .= 2 .* mf.fraction
        # gmf.fraction[mf.γ.interior] .= 2 .* mf.fraction[mf.γ.interior]
        # gmf.fraction[.!mf.γ.interior] .= NaN
        gmf
    end
    analytic_grad = vec(gsoftmax_massfractions(gs, sm))
    # gradients on boundary cells are not touched by gsoftmax_massfractions; treat them as zero
    analytic_grad = map(x -> isnan(x) ? zero(x) : x, analytic_grad)

    println("massfraction softmax cost: ", cost_fn(mvec))
    println("FiniteDiff gradient (first 5): ", fd_grad[1:min(5, length(fd_grad))])
    println("analytic gradient (first 5): ", analytic_grad[1:min(5, length(analytic_grad))])
    @test fd_grad ≈ analytic_grad atol = 1e-8 rtol = 1e-8
    return nothing
end

check_massfraction_softmax_gradient()

@testset "inverse softmax roundtrips" begin
    Random.seed!(42)

    # Vector softmax inverse: re-softmaxing log(s) should give s when s is normalized.
    s = rand(5)
    s ./= sum(s)
    @test softmax_vector(invsoftmax_vector(s)) ≈ s atol = 1e-12 rtol = 1e-12

    # MassFraction softmax inverse: log the softmaxed fractions and recover on wet points.
    sm = softmax_massfractions(m0)
    pre = invsoftmax_massfractions(sm)
    sm_roundtrip = softmax_massfractions(pre)
    for k in keys(sm)
        mask = sm[k].γ.wet .& .!isnan.(sm[k].fraction)
        @test sm_roundtrip[k].fraction[mask] ≈ sm[k].fraction[mask] atol = 1e-12 rtol = 1e-12
    end
end

# Convenience printer to show inverse-softmax roundtrips.
function print_inverse_softmax_roundtrip()
    Random.seed!(123)

    v = randn(5)
    s_from_v = softmax_vector(v)
    inv_from_s_from_v = invsoftmax_vector(s_from_v)
    s_roundtrip = softmax_vector(inv_from_s_from_v)
    println("vector softmax roundtrip:")
    println("  v                         = ", v)
    println("  softmax_vector(v)         = ", s_from_v)
    println("  invsoftmax_vector(softmax(v)) = ", inv_from_s_from_v)
    println("  softmax(inv(softmax(v)))      = ", s_roundtrip)

    sm = softmax_massfractions(m0)
    pre = invsoftmax_massfractions(sm)
    sm_roundtrip = softmax_massfractions(pre)
    k = first(keys(sm))
    I = findfirst(sm[k].γ.wet)
    println("\nmassfraction softmax roundtrip (first wet cell, key $(k)):")
    println("  sm[$(k)].fraction[I]         = ", sm[k].fraction[I])
    println("  invsoftmax_massfractions(sm) = ", pre[k].fraction[I])
    println("  softmax_massfractions(pre)   = ", sm_roundtrip[k].fraction[I])
end

print_inverse_softmax_roundtrip()
