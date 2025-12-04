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
