import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using FiniteDiff

ngrid = (10) # number of grid cells
xmax = 10.0 # domain size 
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

loc = 2.5
δ = interpweights(loc, γ)
@show sum(δ)               # ≈ 1.0
@show sum(δ .* γ.axes[1])       # ≈ loc

sum(δ .* γ.axes)

δ
loc = [2.5]
δ = interpweights(loc, γ)
γ

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


function unconstrained_global_costfunction(control_vector::Vector, 
                                           controls::ControlParameters, c_obs, γ; 
                                           locs = nothing, return_gradients = false)
    du, dq, m = unvec(controls, control_vector)
    if !return_gradients
        J = unconstrained_global_costfunction(du, dq, m, controls, c_obs, γ; 
                                              locs = locs, return_gradients = return_gradients)
        return J
    else
        J, gdu, gdq, gm = unconstrained_global_costfunction(du, dq, m, controls, c_obs, γ; 
                                              locs = locs, return_gradients = return_gradients)
        gcontrols = vcat(vec.([gdu, gdq, gm])...)

        return J, gcontrols
    end
    
end

function unconstrained_global_costfunction(du, dq, m, controls::ControlParameters, 
                                           c_obs, γ; locs = nothing, return_gradients = false)
    # Unpack control vector
    # du, dq, m = unvec(controls, control_vector)
    du_names = keys(du)
    dq_names = keys(dq)

    # Forward pass: compute state at current control vector
    b = deepcopy(controls.u₀)
    q = deepcopy(controls.q₀)

    for key in du_names
        adjustboundarycondition!(b[key], du[key])
    end

    for key in dq_names
        if !isnothing(q[key])
            adjustsource!(q[key], dq[key])
        end
    end

    A = watermassmatrix(m, γ)
    Alu = lu(A)

    # Solve for each tracer
    c = steadyinversion(Alu,b, q, γ)
    n = model_data_misfit(c, c_obs, γ; locs=locs)
    J = model_observation_cost(n,c_obs) + 
        prior_source_cost(dq, controls.q₀, controls.Qₛ) + 
        prior_boundary_cost(du, controls.u₀, controls.Qᵤ) + 
        prior_mass_fraction_cost(m,controls.m₀,controls.Qₘ)
        
    if !return_gradients 
        return J #just return cost
    else
        # Propagate gradients through adjust operations
        gdu = zero(du)
        gdq = zero(dq)
        gm = similar(m0)

        gdu_1 = gprior_boundary_cost(du,controls.u₀, controls.Qᵤ)
        gdq_1 = gprior_source_cost(dq,controls.q₀, controls.Qₛ)
        gm_1 = gprior_mass_fraction_cost(m, controls.m₀, controls.Qₘ)

        # # Adjoint pass: gradient seed is all ones (gradient of sum)
        gn = gmodel_observation_cost(n, c_obs)
        gc = gmodel_data_misfit(gn, c, c_obs, γ; locs=locs)
        # # Compute adjoint gradients
        gdu_2, gdq_2, gA_total = gsteadyinversion(gc, c, A, Alu, b, q, γ)
        gm_2 = gwatermassmatrix(gA_total, m, γ)

        for key in du_names
            gadjustboundarycondition!(gdu[key], gdu_1[key])
            gadjustboundarycondition!(gdu[key], gdu_2[key])
        end

        for key in dq_names
            if !isnothing(q[key])
                gadjustsource!(gdq[key], gdq_1[key], controls.q₀[key])
                gadjustsource!(gdq[key], gdq_2[key], controls.q₀[key])
            end
        end

        for key in keys(gm)
            gm[key].fraction .= gm_1[key].fraction .+ gm_2[key].fraction
        end

        return J, gdu, gdq, gm
    end
end


#deepcopy is necessary to avoid overalpping here 
u₀ = (c = deepcopy(b), c_q = deepcopy(b))
q₀ = (c = nothing, c_q = deepcopy(q),)
c0 = steadyinversion(Alu,u₀, q₀, γ)
Alu = lu(A)

W = (; (k => Diagonal(one.(vec(v))) for (k,v) in pairs(c0))...)
c_obs = (; (k => Observations(1 * v, W[k]) for (k,v) in pairs(c0))...)

du = zero(u₀)
dq = zero(q₀)

diag_from_v(v) = isnothing(v) ? nothing : Diagonal(one.(vec(v)))
Qᵤ = (; (k => diag_from_v(v) for (k,v) in pairs(du))...)
Qₛ = (; (k => diag_from_v(v) for (k,v) in pairs(dq))...)
Qₘ = Diagonal(one.(vec(m0)))

controls = ControlParameters(; du = du, dq = dq, m = m0,
                                u₀ = u₀, q₀ = q₀, m₀ = m0, 
                                Qᵤ = Qᵤ, Qₛ = Qₛ, Qₘ = Qₘ)

                                model_data_misfit

objective(x) = unconstrained_global_costfunction(x, controls, c_obs, γ)

control_vector = randn(length(vec(controls)))
x = randn(length(control_vector))

Jg_finite = FiniteDiff.finite_difference_gradient(objective, control_vector, Val{:central})
J, Jg = unconstrained_global_costfunction(control_vector, controls, c_obs, γ; return_gradients = true)

@. percent_difference(x, y) = 100 * (x - y) / y
all(abs.(percent_difference(Jg_finite, Jg)) .< 0.1) #all within 1 percent 

# using Plots
# scatter(1:length(Jg), Jg)
# scatter!(1:length(Jg), Jg_finite)


