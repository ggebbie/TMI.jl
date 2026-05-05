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
    ngrid = (20,)
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

    m = massfractions_isotropic(γ)
    c = Field(tracer, γ, :c, "linear equilibrated tracer", "μmol/kg")

    b_template = TMI.getboundarycondition(c, 1, 1, γ)
    bvec0 = vec(b_template)

    q_tracer = collect(range(0.0, 0.1, length = ngrid[1]))
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
        grad_m_ad = Enzyme.make_zero(m)
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Reverse),
            J,
            Const(b_template),
            Duplicated(bvec0, grad_b_ad),
            Const(q_template),
            Duplicated(qvec0, grad_q_ad),
            Duplicated(m, grad_m_ad),
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
        grad_m_fd = Enzyme.make_zero(m)

        for k in eachindex(m)
            m1 = m[k]
            for I in TMI.cartesianindex(m1.γ.wet)
                m_plus = deepcopy(m); m_plus[k].fraction[I] += δ
                m_minus = deepcopy(m); m_minus[k].fraction[I] -= δ

                grad_m_fd[k].fraction[I] =
                    (J(b_template, bvec0, q_template, qvec0, m_plus, γ) -
                     J(b_template, bvec0, q_template, qvec0, m_minus, γ)) / (2δ)
            end
        end

        @test all(isfinite.(grad_b_ad))
        @test all(isfinite.(grad_q_ad))
        @test isapprox(grad_b_ad, grad_b_fd; rtol = 1e-4, atol = 1e-8)
        @test isapprox(grad_q_ad, grad_q_fd; rtol = 1e-4, atol = 1e-8)

        for (gad, gfd, m1) in zip(grad_m_ad, grad_m_fd, m)
            mask = m1.γ.wet
            @test all(isfinite.(gad.fraction[mask]))
            @test isapprox(gad.fraction[mask], gfd.fraction[mask]; rtol = 1e-4, atol = 1e-8)
        end
    end
end
