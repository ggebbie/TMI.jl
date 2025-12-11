function centered_finite_difference(f, x; δ = 1e-3)
    orig_size = size(x)
    x_vec = vec(copy(x))
    grad_vec = similar(x_vec)
    for i in eachindex(x_vec)
        x_plus = copy(x_vec)
        x_minus = copy(x_vec)
        x_plus[i] += δ
        x_minus[i] -= δ
        grad_vec[i] = (f(x_plus) - f(x_minus)) / (2δ)
    end
    return reshape(grad_vec, orig_size)
end

percent_difference(a, b) = @. 100 * ((a - b) / b)

@testset "TMI Inversion Tests" begin
    ## Common Setup
    Random.seed!(42)
    
    ngrid = (50,)
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

    m_iso = massfractions_isotropic(γ)
    m0 = (west = m_iso[1], east = m_iso[2])
    
    c_template = Field(tracer, γ, :c, "linear equilibrated tracer", "μmol/kg")
    
    dim = 1
    b_reference = (west = TMI.getboundarycondition(c_template, dim, 1, γ),
                   east = TMI.getboundarycondition(c_template, 1, ngrid[dim], γ))

    qfield = 1.0e-2 * ones(ngrid)

    ## Generic Gradient Check Function Helper
    # This function compares the analytical gradient from the adjoint model
    # with a numerical gradient from finite differences.
    function gradient_check(controls, obs, γ; seed = 1)
        Random.seed!(seed)
        control_vector = randn(length(vec(controls)))
        objective(x) = optim_fg_unconstrained_global_costfunction!(NaN, nothing, x, controls, obs, γ)
        
        fd_grad = centered_finite_difference(objective, control_vector; δ = 1e-4)
        analytic = zero(control_vector)
        optim_fg_unconstrained_global_costfunction!(NaN, analytic, control_vector, controls, obs, γ)
        
        abspdiff = abs.(percent_difference(fd_grad, analytic))
        if !all(abspdiff .< 0.1)
            println(fd_grad[abspdiff .> 0.1])
            println(analytic[abspdiff .> 0.1])
            println(findall(analytic[abspdiff .> 0.1]))
        end
        @test all(abspdiff .< 0.1)
    end

    @testset "inversion gradient checks" begin
        # ### Forward Model and Basic Fields Test
        A = watermassmatrix(m0, γ)
        Alu = lu(A)
        
        q = TMI.Source(-qfield, γ, :q, "remineralized stuff", "μmol/kg", false)

        @testset "steady solutions" begin
            # Test that a steady-state solution for a conservative tracer is correct
            c̃ = steadyinversion(A, b_reference, γ)
            @test maximum(abs.(c̃ - c_template)) < 1e-4

            # Test that a non-conservative tracer with a sink is always lower concentration
            c_noncons = steadyinversion(A, b_reference, γ; q = q)
            Δc = c_noncons - c_template
            @test iszero(sum(Δc .< 0.0))
        end

        # ### Multi-Tracer Inversion Setup
        # Test two tracers:
        # 1. `tracer_boundary`: controlled only by boundary conditions.
        # 2. `tracer_source`: controlled by boundary conditions and an interior source.
        u₀_template = (tracer_boundary = deepcopy(b_reference), tracer_source = deepcopy(b_reference))
        q₀_template = (tracer_boundary = nothing, tracer_source = deepcopy(q))

        # Get a "true" state to test against
        c0 = steadyinversion(Alu, u₀_template, q₀_template, γ)
        Alu = lu(A) # refresh LU factorization after it's been used

        # For the inversion, we will use a zero-value prior for the boundary condition
        u₀_prior = zero(u₀_template)

        boundary_controls = BoundaryControls(u₀_prior;
            variance = (tracer_boundary=1.0, tracer_source=1.0)
        )
        source_controls = SourceControls(q₀_template;
            variance = (tracer_source = 1.0,)
        )
        massfrac_controls = MassFracControls(m0;
            variance = (west=1.0, east=1.0),
            γ = γ
        )
        controls = Controls(γ,
            boundary = boundary_controls,
            source = source_controls,
            massfrac = massfrac_controls
        )

        ## Gradient Check for Observation Types
        @testset "point observation gradient" begin
            obs_loc = [lon[2]]
            observed_vals = map(v -> observe(v, obs_loc, γ), c0)
            W_pt = map(vals -> Diagonal(one.(vals)), observed_vals)
            c_obs = map((vals, w) -> Observations(vals; locs = obs_loc, γ = γ, W = w), observed_vals, W_pt)
            gradient_check(controls, c_obs, γ; seed = 2)
        end

        @testset "full-field observation gradient" begin
            W_full = map(v -> Diagonal(one.(vec(v))), c0)
            c_obs = map((v, w) -> Observations(v; W = w), c0, W_full)
            gradient_check(controls, c_obs, γ; seed = 3)
        end
    end

@testset "source coupling gradient checks" begin
        # This test verifies the adjoint gradient calculation when one tracer's source
        # is a linear function of an independent tracer's source.
        # Here, we model :o2 as dependent on :p1, while :p2 is another independent tracer.
        A = watermassmatrix(m0, γ)
        Alu = lu(A)
        
        source_p1_independent = TMI.Source(-qfield, γ, :p1, "independent p1 source", "μmol/kg", false)
        source_p2_independent = TMI.Source(-qfield .* 2, γ, :p2, "independent p2 source", "μmol/kg", false)
        source_o2_prior_dependent = TMI.Source(-qfield .* 10, γ, :o2, "o2 prior", "μmol/kg", false) # Prior for the dependent tracer
        
        # Priors for all sources that exist in the model
        source_priors_coupled_all = (p1 = deepcopy(source_p1_independent), o2 = deepcopy(source_o2_prior_dependent), p2 = deepcopy(source_p2_independent))
        
        # Priors for boundary conditions, which must include all tracers
        boundary_priors_coupled_all = (p1 = deepcopy(b_reference), o2 = deepcopy(b_reference), p2 = deepcopy(b_reference))
        u₀_prior = zero(boundary_priors_coupled_all) # Use a zero prior for the boundary cost term

        # ### 2. Define initial guess and coupling relationship
        # The initial guess for the controls (`uq`) has `nothing` for the dependent source.
        source_controls_initial_guess_coupled = (p1 = deepcopy(source_p1_independent), o2 = nothing, p2 = deepcopy(source_p2_independent))
        couplings = (o2 = (:p1, -170.0),) # Defines that Source(:o2) = -170.0 * Source(:p1)

        # ### 3. Set up controls using the new helper function
        boundary_controls_coupled = BoundaryControls(u₀_prior;
            ub = boundary_priors_coupled_all,
            variance = (p1=1.0, o2=1.0, p2=1.0)
        )
        source_controls_coupled = SourceControls(source_priors_coupled_all;
            uq = source_controls_initial_guess_coupled,
            dependencies = couplings,
            variance = (p1=1.0, p2=1.0)
        )
        massfrac_controls_coupled = MassFracControls(m0;
            variance = (west=1.0, east=1.0),
            γ = γ
        )
        controls_coupled = Controls(γ;
            boundary = boundary_controls_coupled,
            source = source_controls_coupled,
            massfrac = massfrac_controls_coupled
        )
        
        # Test that the structure of the control parameter gradients is correct for coupled sources
        if !isnothing(couplings) && !isempty(couplings)
            @test length(vec(controls_coupled.source.gduq)) < length(vec(controls_coupled.source.gduq_cache))
        end

        # ### 4. Create synthetic observations for all tracers
        c_p1_true = steadyinversion(Alu, b_reference, γ, q=source_p1_independent)
        c_o2_true = steadyinversion(Alu, b_reference, γ, q=-170.0 * source_p1_independent)
        c_p2_true = steadyinversion(Alu, b_reference, γ, q=source_p2_independent)
        c0_coupled_all = (p1=c_p1_true, o2=c_o2_true, p2=c_p2_true)
        
        # Test that the multi-tracer solve gives the same result as separate solves
        q_all_sources = (p1 = source_p1_independent, o2 = -170.0 * source_p1_independent, p2 = source_p2_independent)
        c_all_at_once = steadyinversion(Alu, boundary_priors_coupled_all, q_all_sources, γ)
        @test isapprox(c0_coupled_all.p1.tracer, c_all_at_once.p1.tracer)
        @test isapprox(c0_coupled_all.o2.tracer, c_all_at_once.o2.tracer)
        @test isapprox(c0_coupled_all.p2.tracer, c_all_at_once.p2.tracer)

        W_coupled = map(v -> Diagonal(one.(vec(v))), c0_coupled_all)
        c_obs_coupled = map((v, w) -> Observations(v; W = w), c0_coupled_all, W_coupled)

        @testset "full-field observation gradient with coupling" begin
            gradient_check(controls_coupled, c_obs_coupled, γ; seed = 4)
        end
    end

    @testset "decay rate feature gradient check" begin
        # This test verifies the adjoint gradient calculation for a tracer that includes
        # a linear decay term. This is useful for modeling radioactive tracers where the
        # decay is proportional to the concentration (λc). The decay is introduced
        # via the `Observations` object for the specific tracer.
        
        A_orig = watermassmatrix(m0, γ)
        Alu_orig = lu(A_orig)

        ### 1. Define priors and initial guesses
        source_p1_independent = TMI.Source(-qfield, γ, :p1, "independent p1 source", "μmol/kg", false)
        source_p1_prior = (p1 = deepcopy(source_p1_independent),)
        
        boundary_p1_initial_guess = (p1 = deepcopy(b_reference),)
        boundary_p1_prior = zero(boundary_p1_initial_guess)

        # ### 2. Set up controls using the new helper function
        boundary_controls_single = BoundaryControls(boundary_p1_prior;
            ub = boundary_p1_initial_guess,
            variance = (p1=1.0,)
        )
        source_controls_single = TMI.SourceControls(source_p1_prior;
            variance = (p1=1.0,)
        )
        massfrac_controls_single = MassFracControls(m0;
            variance = (west=1.0, east=1.0),
            γ = γ
        )
        controls_single_tracer = Controls(γ;
            boundary = boundary_controls_single,
            source = source_controls_single,
            massfrac = massfrac_controls_single
        )

        # ### 3. Create synthetic observations with a decay rate
        c_p1_true = steadyinversion(Alu_orig, b_reference, γ, q=source_p1_independent)
        c0_single = (p1=c_p1_true,)
        W_single = map(v -> Diagonal(one.(vec(v))), c0_single)

        # Define a decay rate and include it in the Observations object.
        # This will modify the transport matrix `A` for this tracer during the solve.
        decay_rate_val = 0.05
        c_obs_decayed_p1 = Observations(c0_single.p1; W = W_single.p1, decay_rate = decay_rate_val, γ=γ)
        c_obs_with_decay = (p1 = c_obs_decayed_p1,)

        # ### 4. Run the gradient check
        gradient_check(controls_single_tracer, c_obs_with_decay, γ; seed = 6)
    end
end