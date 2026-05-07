using Enzyme

function loss(
    b_template::NamedTuple{names, <:Tuple{Vararg{<:BoundaryCondition}}},
    bvec::AbstractVector{T},
    q_template::NamedTuple{qnames, <:Tuple{Vararg{<:TMI.Source}}},
    qvec::AbstractVector{T},
    m::NamedTuple,
    γ::Grid{T},
    b0::NamedTuple{names, <:Tuple{Vararg{<:BoundaryCondition}}},
    q0::NamedTuple{qnames},
    c0::NamedTuple{names, <:Tuple{Vararg{<:Field}}},
) where {names, qnames, T <: Real}
    A = watermassmatrix(m, γ)
    b = unvec(b_template, bvec)
    q = unvec(q_template, qvec)
    ĉ = steadyinversion(A, b, q, γ)

    Δb = b - b0
    Δc = ĉ - c0
    Δq = q - q0

    nΔb = vec(Δb)' * I * vec(Δb)
    nΔc = vec(Δc)' * I * vec(Δc)
    nΔq = vec(Δq)' * I * vec(Δq)
    return nΔb + nΔc + nΔq
end

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
    tracer = fill(1.0, ngrid)
    c = Field(tracer, γ, :c, "linear equilibrated tracer", "umol/kg")

    b_surface = getsurfaceboundary(c)
    b_surface2 = 1.01 * getsurfaceboundary(c)
    b_template = (
        tracer = b_surface,
        tracer2 = b_surface2,
    )
    bvec0 = vec(b_template)

    q = onesource(γ, :q, "interior source", "umol/kg")
    q.tracer[q.γ.interior] .*= 0.01
    q2 = onesource(γ, :q2, "interior source 2", "umol/kg")
    q2.tracer[q2.γ.interior] .*= 0.008
    q_template = (tracer = q, tracer2 = q2)
    qvec0 = vec(q_template)
    q0 = map(qi -> 0.97 * qi, q_template)
    b0 = map(bi -> 0.98 * bi, b_template)
    c0 = map(ci -> 0.95 * ci, steadyinversion(watermassmatrix(m, γ), b0, q0, γ))

    @testset "reverse NamedTuple steadyinversion through unvec and watermassmatrix (b, q, m)" begin
        grad_b_ad = Enzyme.make_zero(bvec0)
        grad_q_ad = Enzyme.make_zero(qvec0)
        grad_m_ad = Enzyme.make_zero(m)
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Reverse),
            loss,
            Const(b_template),
            Duplicated(bvec0, grad_b_ad),
            Const(q_template),
            Duplicated(qvec0, grad_q_ad),
            Duplicated(m, grad_m_ad),
            Const(γ),
            Const(b0),
            Const(q0),
            Const(c0),
        )

        δ = 1e-6
        grad_b_fd = finite_difference_gradient(
            bv -> loss(b_template, bv, q_template, qvec0, m, γ, b0, q0, c0),
            bvec0;
            δ,
        )
        grad_q_fd = finite_difference_gradient(
            qv -> loss(b_template, bvec0, q_template, qv, m, γ, b0, q0, c0),
            qvec0;
            δ,
        )

        m_locs = [(k, I) for k in eachindex(m) for I in TMI.cartesianindex(m[k].γ.wet)]
        mvec0 = [m[k].fraction[I] for (k, I) in m_locs]

        m_from_vec = function (mvec)
            mout = deepcopy(m)
            for (ii, (k, I)) in enumerate(m_locs)
                mout[k].fraction[I] = mvec[ii]
            end
            return mout
        end

        grad_m_fd_vec = finite_difference_gradient(
            mv -> loss(b_template, bvec0, q_template, qvec0, m_from_vec(mv), γ, b0, q0, c0),
            mvec0;
            δ,
        )

        grad_m_fd = Enzyme.make_zero(m)

        for (ii, (k, I)) in enumerate(m_locs)
            grad_m_fd[k].fraction[I] = grad_m_fd_vec[ii]
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
