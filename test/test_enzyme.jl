using Enzyme

"""
    centered_finite_difference(f, x; δ = 1e-6)

Centered FD reference gradient: `dfᵢ ≈ (f(x + δeᵢ) - f(x - δeᵢ)) / 2δ`.
Used to verify the Enzyme gradient.
"""
function centered_finite_difference(f, x; δ = 1e-6)
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
    # ngrid kept small for fast prototyping. Reverse-mode Enzyme is slower than
    # centered FD below ~1000 grid points and faster above (~1.9× at 10_000).
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
        # Scalar objective J = ‖c‖² over wet cells, where:
        #
        #     A = watermassmatrix(m, γ)
        #     b = unvec(b_template, bvec)
        #     q = unvec(q_template, qvec)
        #     c = steadyinversion(A, b, γ; q = q)
        #
        # This single test checks that reverse-mode Enzyme can propagate through
        # both structural reconstruction (`unvec`) and matrix construction
        # (`watermassmatrix`) before the steady solve.
        function loss(b_template, bvec, q_template, qvec, m, γ)
            A = watermassmatrix(m, γ)
            b = unvec(b_template, bvec)
            q = unvec(q_template, qvec)
            ĉ = steadyinversion(A, b, γ; q = q)
            return sum(abs2, ĉ.tracer[ĉ.γ.wet])
        end

        grad_b_ad = Enzyme.make_zero(bvec0)
        grad_q_ad = Enzyme.make_zero(qvec0)
        grad_m_ad = Enzyme.make_zero(m)
        enzyme_time = @elapsed Enzyme.autodiff(
            Enzyme.set_runtime_activity(Reverse),
            loss,
            Const(b_template),
            Duplicated(bvec0, grad_b_ad),
            Const(q_template),
            Duplicated(qvec0, grad_q_ad),
            Duplicated(m, grad_m_ad),
            Const(γ),
        )

        grad_b_fd = similar(bvec0)
        grad_q_fd = similar(qvec0)
        grad_m_fd = Enzyme.make_zero(m)
        finite_difference_time = @elapsed begin
            grad_b_fd = centered_finite_difference(
                bv -> loss(b_template, bv, q_template, qvec0, m, γ),
                bvec0,
            )
            grad_q_fd = centered_finite_difference(
                qv -> loss(b_template, bvec0, q_template, qv, m, γ),
                qvec0,
            )

            # FD reference for `m`: perturb each wet entry of each
            # mass-fraction direction one at a time. Non-wet entries do not
            # enter `A` and have zero gradient by construction, so we skip them.
            δ = 1e-6
            for k in eachindex(m)
                m1 = m[k]
                for I in TMI.cartesianindex(m1.γ.wet)
                    m_plus  = deepcopy(m); m_plus[k].fraction[I]  += δ
                    m_minus = deepcopy(m); m_minus[k].fraction[I] -= δ
                    grad_m_fd[k].fraction[I] =
                        (loss(b_template, bvec0, q_template, qvec0, m_plus, γ) -
                         loss(b_template, bvec0, q_template, qvec0, m_minus, γ)) / (2δ)
                end
            end
        end

        @info "Enzyme reverse-mode timing (b, q, m)" ngrid enzyme_time finite_difference_time speedup = finite_difference_time / enzyme_time
        @test all(isfinite.(grad_b_ad))
        @test all(isfinite.(grad_q_ad))
        # isapprox handles entries where the FD reference is ~0 (percent
        # difference would divide by zero there).
        @test isapprox(grad_b_ad, grad_b_fd; rtol = 1e-4, atol = 1e-8)
        @test isapprox(grad_q_ad, grad_q_fd; rtol = 1e-4, atol = 1e-8)
        for (gad, gfd, m1) in zip(grad_m_ad, grad_m_fd, m)
            mask = m1.γ.wet
            @test all(isfinite.(gad.fraction[mask]))
            @test isapprox(gad.fraction[mask], gfd.fraction[mask]; rtol = 1e-4, atol = 1e-8)
        end
    end
end
