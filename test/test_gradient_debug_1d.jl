
import Pkg; Pkg.activate(".")

using Test
using TMI
using LinearAlgebra
using Optim
using FiniteDiff # Changed from ForwardDiff
using LinearSolve, IncompleteLU

@testset "Gradient Debugging with FiniteDiff" begin # Changed testset name

    # 1. Set up a simple 1D problem
    ngrid = (10,)
    xmax = 1000.0
    lon = collect(range(0.0, xmax, length=ngrid[1]))
    
    axes = (lon,)
    wet = trues(ngrid)
    interior = copy(wet)
    interior[begin] = false
    interior[end] = false
    wrap = (false,)
    Δ = [CartesianIndex(1,), CartesianIndex(-1,)]
    γ = Grid(axes, wet, interior, wrap, Δ)

    # 2. Create observations and priors for two tracers
    c_field = Field(collect(1.0 .- lon./xmax), γ, :c, "linear tracer", "μmol/kg")
    s_field = Field(collect((lon./xmax).^2), γ, :S, "quadratic tracer", "psu")
    cobs = (c = c_field, S = s_field)
    
    dim = 1
    # Boundary conditions for the `c` tracer at the west and east ends
    b_c_tracer = (west = TMI.getboundarycondition(c_field, dim, 1, γ),
                  east = TMI.getboundarycondition(c_field, dim, ngrid[dim], γ))
    b_s_tracer = (west = TMI.getboundarycondition(s_field, dim, 1, γ),
                  east = TMI.getboundarycondition(s_field, dim, ngrid[dim], γ))

    u₀ = (c = b_c_tracer, S = b_s_tracer)
    q₀ = map(v -> nothing, cobs)
    
    m_iso = massfractions_isotropic(γ)
    m₀ = (west = m_iso[1], east = m_iso[2])

    # Manual setup of controls and priors, similar to total_inversion_3d.jl
    ub = deepcopy(u₀)
    uq = deepcopy(q₀)
    m = deepcopy(m₀)

    # Variances
    w = (c = 0.01^2, S = 0.1^2)
    m_variance = (west = 1.0, east = 1.0)

    # create a vector for mass fraction variances
    m_variance_vec = vcat([fill(m_variance[k], length(vec(m₀[k]))) for k in keys(m₀)]...)

    # Precision matrices
    Qᵤ = map((v, wi) -> Diagonal(one.(vec(v)) .* 1/wi), ub, w)
    Qₛ = map(v -> nothing, uq)
    Qₘ = Diagonal( 1.0 ./ m_variance_vec )

    controls = ControlParameters(; γ = γ, ub = ub, uq = uq, m = m,
                                u₀ = u₀, q₀ = q₀, m₀ = m₀, 
                                Qᵤ = Qᵤ, Qₛ = Qₛ, Qₘ = Qₘ)

    # Use unit weights for observations, variance is handled in the prior cost
    W_full = map(v -> Diagonal(one.(vec(v))), cobs)
    c_obs_ad = map((v, w) -> Observations(v; W = w), cobs, W_full)

    # 3. Define the cost function for FiniteDiff
    # This function must not be in-place
    function cost_function_finitediff(x_vec::Vector{T}) where T # Changed function name
        
        # Create a new controls object to avoid modifying the global one
        local_controls = deepcopy(controls)

        # Unvectorize the control vector
        unvec!(local_controls, x_vec)

        # Apply softmax
        α = 5.0
        softmax_massfractions!(local_controls.massfrac.m; α=α)

        # Run the forward pass
        state = TMI.unconstrained_global_forward(local_controls, c_obs_ad, γ)
        
        return state.J
    end

    # 4. Create an initial control vector
    # Use invsoftmax to start near the prior
    m_inv_sm = invsoftmax_massfractions(m₀; α = 5.0)
    x0_m = vec(m_inv_sm)
    
    # Use priors for boundary conditions
    x0_u = vec(u₀)
    
    x0 = vcat(x0_u, x0_m)

    # 5. Calculate gradients
    
    # Gradient from FiniteDiff
    grad_finitediff = FiniteDiff.finite_difference_gradient(cost_function_finitediff, x0) # Changed call

    # Gradient from the analytical function
    grad_analytical = similar(x0)
    fg!(F, G, x) = optim_fg_constrained_global_costfunction!(F, G, x, controls, c_obs_ad, γ)
    fg!(nothing, grad_analytical, x0)

    # 5b. Quick gradient check on a random element
    println("\nPerforming single-element centered finite difference check...")
    
    idx = rand(1:length(x0))
    small = 1e-8
    
    x_plus = deepcopy(x0); x_plus[idx] += small
    x_minus = deepcopy(x0); x_minus[idx] -= small
    
    grad_finite_element = (cost_function_finitediff(x_plus) - cost_function_finitediff(x_minus)) / (2 * small)
    grad_analytical_element = grad_analytical[idx]
    
    println("Random index [", idx, "]:")
    println("  Analytical:  ", grad_analytical_element)
    println("  Finite Diff: ", grad_finite_element)
    @test grad_analytical_element ≈ grad_finite_element atol=1e-6

    # 6. Compare the gradients
    println("\nComparing analytical gradient to FiniteDiff gradient...") # Changed print statement
    println("Norm of difference: ", norm(grad_analytical - grad_finitediff))
    println("Relative norm of difference: ", norm(grad_analytical - grad_finitediff) / norm(grad_analytical))
    
    @test grad_analytical ≈ grad_finitediff atol=1e-8 # Changed grad_forwarddiff to grad_finitediff
    
end
