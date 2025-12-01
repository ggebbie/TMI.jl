#!/usr/bin/env julia
# Benchmark and validate an in-place scratch version of unconstrained_global_costfunction.
# This does NOT modify the TMI source; it defines local helpers and runs checks.

import Pkg; Pkg.activate(".")

using TMI, LinearAlgebra, Statistics, BenchmarkTools, FiniteDiff
using TMI: Source, BoundaryCondition, Field, unconstrained_global_costfunction,
    adjustboundarycondition!, adjustsource!, gsteadyinversion, gwatermassmatrix,
    prior_source_cost, prior_boundary_cost, prior_mass_fraction_cost,
    gprior_boundary_cost, gprior_source_cost, gprior_mass_fraction_cost,
    model_data_misfit, model_observation_cost, gmodel_observation_cost, gmodel_data_misfit
# Local helper to reset a Source from another Source (prototype only).
setsource!(d::Source{T}, q::Source{T}) where T <: Real = (d.tracer[d.γ.interior] .= q.tracer[q.γ.interior]; d)

function zero_source!(s::Source)
    s.tracer[s.γ.interior] .= 0
    return s
end

function make_unconstrained_scratch(controls::ControlParameters, γ)
    b  = deepcopy(controls.u₀)
    q  = deepcopy(controls.q₀)
    du = deepcopy(controls.u)
    dq = deepcopy(controls.q)
    gdu = zero(deepcopy(controls.u))
    gdq = zero(deepcopy(controls.q))
    gm  = map(m -> (mf = deepcopy(m); mf.fraction .= 0; mf), controls.m)
    return (; b, q, du, dq, gdu, gdq, gm)
end

function unconstrained_global_costfunction!(u, q, m, controls::ControlParameters, c_obs, γ;
    scratch,
    locs = nothing,
    return_gradients::Bool = false)

    # reset control perturbations to current ub/uq
    for key in keys(scratch.du)
        setboundarycondition!(scratch.du[key], u[key])
    end
    for key in keys(scratch.dq)
        if !isnothing(scratch.dq[key])
            setsource!(scratch.dq[key], q[key])
        end
    end

    # reset priors into b/q
    for key in keys(scratch.b)
        setboundarycondition!(scratch.b[key], controls.u₀[key])
    end
    for key in keys(scratch.q)
        if !isnothing(scratch.q[key])
            setsource!(scratch.q[key], controls.q₀[key])
        end
    end

    # zero gradients
    for key in keys(scratch.gdu)
        setboundarycondition!(scratch.gdu[key], zero(scratch.gdu[key]))
    end
    for key in keys(scratch.gdq)
        if !isnothing(scratch.gdq[key])
            zero_source!(scratch.gdq[key])
        end
    end
    for key in keys(scratch.gm)
        scratch.gm[key].fraction .= 0
    end

    # Apply control updates: du = ub - u₀; b = u₀ + du; same for dq/q
    for key in keys(scratch.b)
        if key ∈ keys(u)
            adjustboundarycondition!(scratch.du[key], controls.u₀[key]; r = -1.0)
            adjustboundarycondition!(scratch.b[key], scratch.du[key]; r = 1.0)
        end
    end
    for key in keys(scratch.q)
        if !isnothing(scratch.q[key]) && key ∈ keys(q)
            adjustsource!(scratch.dq[key], controls.q₀[key]; r = -1.0)
            adjustsource!(scratch.q[key], scratch.dq[key]; r = 1.0)
        end
    end

    A = watermassmatrix( m, γ)
    Alu = lu(A)
    c = steadyinversion(Alu, scratch.b, scratch.q, γ)
    n = model_data_misfit(c, c_obs, γ; locs = locs)
    J = model_observation_cost(n, c_obs) +
        prior_source_cost(scratch.dq, controls.q₀, controls.Qₛ) +
        prior_boundary_cost(scratch.du, controls.u₀, controls.Qᵤ) +
        prior_mass_fraction_cost(m, controls.m₀, controls.Qₘ)

    if !return_gradients
        return J, nothing, nothing, nothing
    end

    gdu = gprior_boundary_cost(scratch.du, controls.u₀, controls.Qᵤ)
    gdq = gprior_source_cost(scratch.dq, controls.q₀, controls.Qₛ)
    gm = gprior_mass_fraction_cost(m, controls.m₀, controls.Qₘ)
    for key in keys(scratch.gdu)
        setboundarycondition!(scratch.gdu[key], gdu[key])
    end
    for key in keys(scratch.gdq)
        if !isnothing(scratch.gdq[key])
            setsource!(scratch.gdq[key], gdq[key])
        end
    end
    for key in keys(scratch.gm)
        scratch.gm[key].fraction .= gm[key].fraction
    end

    gn = gmodel_observation_cost(n, c_obs)
    gc = gmodel_data_misfit(gn, c, c_obs, γ; locs = locs)
    gdu, gdq, gA_total = gsteadyinversion(gc, c, A, Alu, scratch.b, scratch.q, γ)
    
    for key in keys(scratch.gdu)
        setboundarycondition!(scratch.gdu[key], gdu[key])
    end
    
    for key in keys(scratch.gdq)
        if !isnothing(scratch.gdq[key])
            setsource!(scratch.gdq[key], gdq[key])
        end
    end
   
    gm_dyn = gwatermassmatrix(gA_total, m, γ)
    for key in keys(scratch.gm)
        scratch.gm[key].fraction .+= gm_dyn[key].fraction
    end

    return J, scratch.gdu, scratch.gdq, scratch.gm
end

# --- test setup like scripts/total_inversion.jl (1D) ---
ngrid = (1000,)
xmax = 1000.0
lon = collect(range(0.0, xmax, length = ngrid[1]))
tracer = collect(1.0 .- lon ./ xmax)
axes = (lon,)
wet = trues(ngrid)
interior = copy(wet)
interior[begin] = false
interior[end] = false
wrap = (false,)
Δ = [CartesianIndex(1,), CartesianIndex(-1,)]
γ = Grid(axes, wet, interior, wrap, Δ)
m0 = massfractions_isotropic(γ)
m0 = (west = m0[1], east = m0[2])
c = Field(tracer, γ, :c, "linear equilibrated tracer", "μmol/kg")
A = watermassmatrix(m0, γ)
Alu = lu(A)
dim = 1
b = (west = TMI.getboundarycondition(c, dim, 1, γ), east = TMI.getboundarycondition(c, 1, ngrid[dim], γ))
qfield = 1.0e-2 * ones(ngrid)
q = TMI.Source(-qfield, γ, :q, "remineralized stuff", "μmol/kg", false)
u₀ = (c = deepcopy(b), c_q = deepcopy(b))
q₀ = (c = nothing, c_q = deepcopy(q))
c0 = steadyinversion(Alu, u₀, q₀, γ)
Alu = lu(A)
ub = deepcopy(u₀)
uq = deepcopy(q₀)
Qᵤ = map(v -> Diagonal(one.(vec(v))), ub)
Qₛ = map(v -> isnothing(v) ? nothing : Diagonal(one.(vec(v))), uq)
Qₘ = Diagonal(one.(vec(m0)))
controls = ControlParameters(; u = ub, q = uq, m = m0, u₀ = u₀, q₀ = q₀, m₀ = m0, Qᵤ = Qᵤ, Qₛ = Qₛ, Qₘ = Qₘ)
W_full = map(v -> Diagonal(one.(vec(v))), c0)
c_obs = map((v, w) -> Observations(v; W = w), c0, W_full)

scratch = make_unconstrained_scratch(controls, γ)

u_test = deepcopy(u₀)
q_test = deepcopy(q₀)
m_test = deepcopy(m0)
control_vector = vcat(vec(u_test), vec(q_test), vec(m_test))
objective(x) = unconstrained_global_costfunction(x, controls, c_obs, γ; return_gradients = false)
Jfd, gradfd = FiniteDiff.finite_difference_gradient(objective, control_vector, Val{:central})
Jb, gb_u, gb_q, gb_m = unconstrained_global_costfunction(u_test, q_test, m_test, controls, c_obs, γ; locs = nothing, return_gradients = true)
Js, gs_u, gs_q, gs_m = unconstrained_global_costfunction!(u_test, q_test, m_test, controls, c_obs, γ; scratch = scratch, locs = nothing, return_gradients = true)
gb_vec = vcat(vec(gb_u), vec(gb_q), vec(gb_m))
gs_vec = vcat(vec(gs_u), vec(gs_q), vec(gs_m))
println("Cost diff baseline vs FD: ", abs(Jb - Jfd))
println("Grad max diff baseline vs FD: ", maximum(abs.(gb_vec .- gradfd)))
println("Cost diff scratch vs baseline: ", abs(Js - Jb))
println("Grad max diff scratch vs baseline: ", maximum(abs.(gs_vec .- gb_vec)))

println("Baseline timing (@btime):")
@btime unconstrained_global_costfunction($u_test, $q_test, $m_test, controls, c_obs, γ; locs = nothing, return_gradients = true);
println("Scratch timing (@btime):")
@btime unconstrained_global_costfunction!($u_test, $q_test, $m_test, controls, c_obs, γ; scratch = $scratch, locs = nothing, return_gradients = true);

println("Baseline timing (@btime):")
@btime unconstrained_global_costfunction($u_test, $q_test, $m_test, controls, c_obs, γ; locs = nothing, return_gradients = true);
println("Scratch timing (@btime):")
@btime unconstrained_global_costfunction!($u_test, $q_test, $m_test, controls, c_obs, γ; scratch = $scratch, locs = nothing, return_gradients = true);
