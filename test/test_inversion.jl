using Printf
percent_difference(a, b) = @. 100 * ((a - b) / b)

function centered_finite_difference(f, x; δ = 1e-6)
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
    reshape(grad_vec, orig_size)
end

function generate_1d_grid()
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
    qfield[γ.wet .& .!γ.interior] .= 0.0
    return (γ, lon, c_template, b_reference, qfield, m0)
end

function gradient_check(controls, obs, γ; seed = 1)
    Random.seed!(seed)
    control_vector = randn(length(vec(controls)))
    objective(x) = joint_global_cost!(NaN, nothing, x, controls, obs, γ)

    fd_grad = centered_finite_difference(objective, control_vector; δ = 1e-3)
    analytic = zero(control_vector)
    joint_global_cost!(NaN, analytic, control_vector, controls, obs, γ)

    abspdiff = abs.(percent_difference(fd_grad, analytic))

    # --- DEBUGGING ADDITION ---
    if !all(abspdiff .< 0.1)
        println("DEBUG: Gradient check failed. Mismatch > 0.1%")
        println("Index | Finite Diff | Analytic | % Diff")
        println("----------------------------------------------")
        for i in eachindex(abspdiff)
            if abspdiff[i] >= 0.1
                println(i)
                @printf("  %3d | %11.6f | %11.6f | %.2f%%\n", i, fd_grad[i], analytic[i], abspdiff[i])
            end
        end
    end
    # --- END DEBUGGING ADDITION ---
    
    @test all(abspdiff .< 0.1)
end

@testset "TMI Inversion Tests" begin
    Random.seed!(42)
    γ, lon, c_template, b_reference, qfield, m0 = generate_1d_grid()

    @testset "inversion gradient checks" begin
        A = watermassmatrix(m0, γ)
        Alu = lu(A)

        q = TMI.Source(-qfield, γ, :q, "remineralized stuff", "μmol/kg", false)

        @testset "steady solutions" begin
            c̃ = steadyinversion(A, b_reference, γ)
            @test maximum(abs.(c̃ - c_template)) < 1e-4

            c_noncons = steadyinversion(A, b_reference, γ; q = q)
            Δc = c_noncons - c_template
            @test iszero(sum(Δc .< 0.0))
        end

        u₀_template = (tracer_boundary = deepcopy(b_reference), tracer_source = deepcopy(b_reference))
        q₀_template = (tracer_boundary = nothing, tracer_source = deepcopy(q))

        c0 = steadyinversion(Alu, u₀_template, q₀_template, γ)
        Alu = lu(A) # refresh LU factorization after it's been used

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
        A = watermassmatrix(m0, γ)
        Alu = lu(A)

        source_p1_independent = TMI.Source(-qfield, γ, :p1, "independent p1 source", "μmol/kg", false)
        source_p2_independent = TMI.Source(-qfield .* 2, γ, :p2, "independent p2 source", "μmol/kg", false)
        source_o2_prior_dependent = TMI.Source(-qfield .* 10, γ, :o2, "o2 prior", "μmol/kg", false)

        source_priors_coupled_all = (p1 = deepcopy(source_p1_independent), o2 = deepcopy(source_o2_prior_dependent), p2 = deepcopy(source_p2_independent))
        boundary_priors_coupled_all = (p1 = deepcopy(b_reference), o2 = deepcopy(b_reference), p2 = deepcopy(b_reference))
        u₀_prior = zero(boundary_priors_coupled_all)

        source_controls_initial_guess_coupled = (p1 = deepcopy(source_p1_independent), o2 = nothing, p2 = deepcopy(source_p2_independent))
        couplings = (o2 = (:p1, -170.0),)

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

        if !isnothing(couplings) && !isempty(couplings)
            @test length(vec(controls_coupled.source.gduq)) < length(vec(controls_coupled.source.gduq_cache))
        end

        c_p1_true = steadyinversion(Alu, b_reference, γ, q=source_p1_independent)
        c_o2_true = steadyinversion(Alu, b_reference, γ, q=-170.0 * source_p1_independent)
        c_p2_true = steadyinversion(Alu, b_reference, γ, q=source_p2_independent)
        c0_coupled_all = (p1=c_p1_true, o2=c_o2_true, p2=c_p2_true)

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
end
