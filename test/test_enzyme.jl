using Enzyme

function finite_difference_gradient(f, x; δ = 1e-6)
    grad = similar(x)
    for i in eachindex(x)
        x_plus = copy(x)
        x_minus = copy(x)
        x_plus[i] += δ
        x_minus[i] -= δ
        grad[i] = (f(x_plus) - f(x_minus)) / (2δ)
    end
    return grad
end

@testset "Enzyme Reverse Readiness" begin
    ngrid = (4, 3, 3)
    lon = collect(range(0.0, 1000.0, length = ngrid[1]))
    lat = collect(range(0.0, 800.0, length = ngrid[2]))
    dep = collect(range(0.0, 500.0, length = ngrid[3]))

    axes = (lon, lat, dep)
    wet = trues(ngrid)
    interior = copy(wet)
    interior[:, :, 1] .= false
    interior[:, :, end] .= false
    interior[:, 1, :] .= false
    interior[:, end, :] .= false
    interior[1, :, :] .= false
    interior[end, :, :] .= false

    wrap = (false, false, false)
    Δ = [
        CartesianIndex(1, 0, 0), CartesianIndex(-1, 0, 0),
        CartesianIndex(0, 1, 0), CartesianIndex(0, -1, 0),
        CartesianIndex(0, 0, 1), CartesianIndex(0, 0, -1),
    ]
    γ = Grid(axes, wet, interior, wrap, Δ)

    m = massfractions_isotropic(γ)
    mvec0 = vec(m)
    tracer = fill(1.0, ngrid)
    c = Field(tracer, γ, :c, "linear equilibrated tracer", "μmol/kg")

    b_template = TMI.getboundarycondition(c, 1, 1, γ)
    bvec0 = vec(b_template)

    q_tracer = fill(0.01, ngrid)
    q_template = Field(q_tracer, γ, :q, "interior source", "μmol/kg")
    qvec0 = vec(q_template)

    @testset "reverse steadyinversion through unvec and watermassmatrix (b, q, m)" begin
        function J(b_template, bvec, q_template, qvec, m, γ)
            A = watermassmatrix(m, γ)
            b = unvec(b_template, bvec)
            q = unvec(q_template, qvec)
            ĉ = steadyinversion(A, b, γ; q = q)
            return sum(abs2, ĉ.tracer[ĉ.γ.wet])
        end

        grad_b_ad = Enzyme.make_zero(bvec0)
        grad_q_ad = Enzyme.make_zero(qvec0)
        grad_m_ad = Enzyme.make_zero(mvec0)
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Reverse),
            J,
            Const(b_template),
            Duplicated(bvec0, grad_b_ad),
            Const(q_template),
            Duplicated(qvec0, grad_q_ad),
            Const(m),
            Duplicated(mvec0, grad_m_ad),
            Const(γ),
        )

        δ = 1e-6
        grad_b_fd = finite_difference_gradient(
            bv -> J(b_template, bv, q_template, qvec0, m, γ),
            bvec0;
            δ,
        )
        grad_q_fd = finite_difference_gradient(
            qv -> J(b_template, bvec0, q_template, qv, m, γ),
            qvec0;
            δ,
        )
        grad_m_fd_vec = finite_difference_gradient(
            mv -> J(b_template, bvec0, q_template, qvec0, m_from_vec(mv), γ),
            mvec0;
            δ,
        )

        @test all(isfinite.(grad_b_ad))
        @test all(isfinite.(grad_q_ad))
        @test all(isfinite.(grad_m_ad))
        @test isapprox(grad_b_ad, grad_b_fd; rtol = 1e-4, atol = 1e-8)
        @test isapprox(grad_q_ad, grad_q_fd; rtol = 1e-4, atol = 1e-8)
        @test isapprox(grad_m_ad, grad_m_fd_vec; rtol = 1e-4, atol = 1e-8)
    end
end
