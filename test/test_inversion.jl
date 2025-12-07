function centered_finite_difference(f, x; δ = 1e-12)
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

percent_difference(a, b; eps = 1e-12) = @. 100 * ((a - b) / b)

@testset "inversion gradient checks" begin
    Random.seed!(42)

    ngrid = (50,) # number of grid cells
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
    c = Field(tracer, γ, :c, "linear equilibrated tracer", "μmol/kg")

    A = watermassmatrix(m0, γ)
    Alu = lu(A)

    dim = 1
    b = (west = TMI.getboundarycondition(c, dim, 1, γ),
         east = TMI.getboundarycondition(c, 1, ngrid[dim], γ))

    qfield = 1.0e-2 * ones(ngrid)
    q = TMI.Source(-qfield, γ, :q, "remineralized stuff", "μmol/kg", false)

    @testset "steady solutions" begin
        c̃ = steadyinversion(A, b, γ)
        @test maximum(abs.(c̃ - c)) < 1e-4

        c_noncons = steadyinversion(A, b, γ; q = q)
        Δc = c_noncons - c
        @test iszero(sum(Δc .< 0.0))
    end

    # Controls and priors -------------------------------------------------------
    u₀ = (c = deepcopy(b), c_q = deepcopy(b))
    q₀ = (c = nothing, c_q = deepcopy(q))
    c0 = steadyinversion(Alu, u₀, q₀, γ)
    Alu = lu(A) # refresh LU to avoid reuse of mutated factorization
    TMI.zero!(u₀)

    ub = deepcopy(u₀)
    uq = deepcopy(q₀)
    Qᵤ = map(v -> Diagonal(one.(vec(v))), ub)
    Qₛ = map(v -> isnothing(v) ? nothing : Diagonal(one.(vec(v))), uq)
    Qₘ = Diagonal(one.(vec(m0)))

    controls = ControlParameters(; γ = γ,
                                 ub = ub, uq = uq, m = m0,
                                 u₀ = u₀, q₀ = q₀, m₀ = m0,
                                 Qᵤ = Qᵤ, Qₛ = Qₛ, Qₘ = Qₘ)

    function gradient_check(obs; seed = 1, atol = 1e-10)
        Random.seed!(seed)
        control_vector = randn(length(vec(controls)))
        objective(x) = optim_fg_constrained_global_costfunction!(NaN, nothing, x, controls, obs, γ)
        fd_grad = centered_finite_difference(objective, control_vector; δ = 1e-6)
        analytic = zero(control_vector)
        optim_fg_constrained_global_costfunction!(NaN, analytic, control_vector, controls, obs, γ)
        abspdiff = abs.(percent_difference(fd_grad, analytic; eps = atol))
        @test all(abspdiff .< 0.1) # within 0.1% per component
    end

    @testset "point observation gradient" begin
        obs_loc = [lon[2]]
        observed_vals = map(v -> observe(v, obs_loc, γ), c0)
        W_pt = map(vals -> Diagonal(one.(vals)), observed_vals)
        c_obs = map((vals, w) -> Observations(vals; locs = obs_loc, γ = γ, W = w), observed_vals, W_pt)
        gradient_check(c_obs; seed = 2)
    end

    @testset "full-field observation gradient" begin
        W_full = map(v -> Diagonal(one.(vec(v))), c0) #W is just 1.0 on diagonals for all tracers
        c_obs = map((v, w) -> Observations(v; W = w), c0, W_full)
        gradient_check(c_obs; seed = 3)
    end
end

@testset "source coupling gradient checks" begin
    Random.seed!(42)

    ngrid = (50,) # number of grid cells
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
    c_template = Field(tracer, γ, :c_template, "linear equilibrated tracer", "μmol/kg")

    A = watermassmatrix(m0, γ)
    Alu = lu(A)

    dim = 1
    b_reference = (west = TMI.getboundarycondition(c_template, dim, 1, γ),
                  east = TMI.getboundarycondition(c_template, 1, ngrid[dim], γ))

    qfield = 1.0e-2 * ones(ngrid)
    
    # --- Setup for Coupled Sources ---
    # Independent source tracer :p1 (e.g., PO4)
    q_p1 = TMI.Source(-qfield, γ, :p1, "independent p1 source", "μmol/kg", false)
    # Another independent source tracer :p2
    q_p2 = TMI.Source(-qfield .* 2, γ, :p2, "independent p2 source", "μmol/kg", false)

    # uq_coupled defines the independent controls. :o2 is `nothing`.
    uq_coupled = (p1 = deepcopy(q_p1), o2 = nothing, p2 = deepcopy(q_p2))
    
    # q₀_coupled defines the priors for all tracers that might have sources
    # Even if :o2 is dependent, it needs a prior to derive its `q` buffer.
    q_o2_prior = TMI.Source(-qfield .* 10, γ, :o2, "o2 prior", "μmol/kg", false)
    q₀_coupled_all = (p1 = deepcopy(q_p1), o2 = deepcopy(q_o2_prior), p2 = deepcopy(q_p2))
    
    # Covariances are only needed for the independent controls
    Qₛ_coupled = (p1 = Diagonal(one.(vec(q_p1))), p2 = Diagonal(one.(vec(q_p2))))

    # Define the coupling relationship: Source(o2) = -170.0 * Source(p1)
    couplings = (o2 = (:p1, -170.0),)

    # --- Other controls ---
    # Boundary controls must have the same tracer keys as the sources
    ub = (p1 = deepcopy(b_reference), o2 = deepcopy(b_reference), p2 = deepcopy(b_reference))
    u₀ = zero(ub) # Zero prior for boundary perturbations
    Qᵤ = map(v -> Diagonal(one.(vec(v))), ub)
    Qₘ = Diagonal(one.(vec(m0)))
    m_full = (west = m_iso[1], east = m_iso[2])

    # 3. INSTANTIATE ControlParameters WITH THE NEW ARGUMENT
    controls_coupled = ControlParameters(; γ = γ,
                                 ub = ub, uq = uq_coupled, m = m_full,
                                 u₀ = u₀, q₀ = q₀_coupled_all, m₀ = m0,
                                 Qᵤ = Qᵤ, Qₛ = Qₛ_coupled, Qₘ = Qₘ,
                                 source_couplings = couplings)

    # TEST: length(gduq) < length(gduq_cache) if couplings exist
    if !isnothing(couplings) && !isempty(couplings)
        @test length(vec(controls_coupled.source.gduq)) < length(vec(controls_coupled.source.gduq_cache))
    end

    # 4. CREATE SYNTHETIC OBSERVATIONS FOR ALL TRACERS
    #    We need observations for p1, o2, p2.
    #    First, solve for a steady state to get some "true" values
    c_p1_true = steadyinversion(Alu, b_reference, γ, q=q_p1)
    c_o2_true = steadyinversion(Alu, b_reference, γ, q=-170.0 * q_p1)
    c_p2_true = steadyinversion(Alu, b_reference, γ, q=q_p2)

    c0_coupled_all = (p1=c_p1_true, o2=c_o2_true, p2=c_p2_true)
    
    # New test: Compare with single multi-tracer solve
    q_all_sources = (p1 = q_p1, o2 = -170.0 * q_p1, p2 = q_p2)
    b_all_tracers = (p1 = b_reference, o2 = b_reference, p2 = b_reference)
    c_all_at_once = steadyinversion(Alu, b_all_tracers, q_all_sources, γ)
    @test isapprox(c0_coupled_all.p1.tracer, c_all_at_once.p1.tracer)
    @test isapprox(c0_coupled_all.o2.tracer, c_all_at_once.o2.tracer)
    @test isapprox(c0_coupled_all.p2.tracer, c_all_at_once.p2.tracer)

    W_coupled = map(v -> Diagonal(one.(vec(v))), c0_coupled_all)
    c_obs_coupled = map((v, w) -> Observations(v; W = w), c0_coupled_all, W_coupled)

    # 5. RUN THE GRADIENT CHECK
    function gradient_check_coupled(obs; seed = 1, atol = 1e-10)
        Random.seed!(seed)
        # The control_vector should ONLY include independent sources (p1, p2)
        control_vector = randn(length(vec(controls_coupled)))
        
        objective(x) = optim_fg_constrained_global_costfunction!(
            NaN, nothing, x, controls_coupled, obs, γ
        )
        fd_grad = centered_finite_difference(objective, control_vector; δ = 1e-6)
        analytic = zero(control_vector)
        optim_fg_constrained_global_costfunction!(
            NaN, analytic, control_vector, controls_coupled, obs, γ
        )
        abspdiff = abs.(percent_difference(fd_grad, analytic; eps = atol))
        @test all(abspdiff .< 0.1) # within 0.1% per component
    end

    @testset "full-field observation gradient with coupling" begin
        gradient_check_coupled(c_obs_coupled; seed = 4)
    end
end
