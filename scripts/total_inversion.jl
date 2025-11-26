import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using FiniteDiff

ngrid = (10) # number of grid cells
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

# Alu \ setsource q
function forward_model(control_vector, controls, γ)

    du, dq, m = unvec(controls, control_vector)
    du_names = keys(du)
    dq_names = keys(dq)

    b = deepcopy(controls.u₀)
    q = deepcopy(controls.q₀)

    #might already exist!
    for key in du_names #only update if in the control vector
        adjustboundarycondition!(b[key],du[key]) #b += u # easy case where u and b are on the same boundary
    end

    for key in dq_names
        if !isnothing(dq[key])
            adjustsource!(q[key],dq[key]) #b += u # easy case where u and b are on the same boundary
        end
    end

    A = watermassmatrix(m, γ)
    Alu = lu(A)

    c = steadyinversion(Alu,b, q, γ)

    return sum(vec(c))
    # return c
end

function adjoint_gradient(control_vector, controls, γ)
    """
    Compute adjoint gradient by first running forward model to get state,
    then running adjoint code.
    """
    # Unpack control vector
    du, dq, m = unvec(controls, control_vector)
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

    # Adjoint pass: gradient seed is all ones (gradient of sum)
    gc =one(c)

    # Compute adjoint gradients
    gb, gq, gA_total = gsteadyinversion(gc, c, A, Alu, b, q, γ)
    gm = gwatermassmatrix(gA_total, m, γ)

    # # Propagate gradients through adjust operations
    # gub = zero(gb)
    # guq = zero(gq)

    # for key in du_names
    #     gadjustboundarycondition!(gub[key], gb[key])
    # end

    # for key in dq_names
    #     if !isnothing(q[key])
    #         gadjustsource!(guq[key], gq[key])
    #     end
    # end

    return gb, gq, gm
end


#deepcopy is necessary to avoid overalpping here 
u₀ = (c = deepcopy(b), c_q = deepcopy(b))
q₀ = (c = nothing, c_q = deepcopy(q),)
c0 = steadyinversion(Alu,u₀, q₀, γ)
Alu = lu(A)


du = zero(u₀)
dq = zero(q₀)
Qⁱᵤ = Diagonal(one.(vec(du)))
Qⁱₛ = Diagonal(one.(vec(dq)))
Qⁱₘ = Diagonal(one.(vec(m0)))


controls = ControlParameters(; du = du, dq = dq, m = m0,
                                u₀ = u₀, q₀ = q₀, m₀ = m0)

objective(x) = forward_model(x, controls, γ)
control_vector = vcat(vec.([du,dq, m0])...)
g = FiniteDiff.finite_difference_gradient(objective, control_vector, Val{:central})
gub, guq, gm = adjoint_gradient(control_vector, controls, γ)
g_adjoint = vectorize_controls(gub, guq, gm) #boundary condition
all(g_adjoint .≈ g)