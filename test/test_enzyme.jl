using Enzyme
using LinearAlgebra: I

function loss(
    bvec::AbstractVector{T},
    qvec::AbstractVector{T},
    mvec::AbstractVector{T},
    γ::Grid{T},
    b0::NamedTuple{names, <:Tuple{Vararg{<:BoundaryCondition}}},
    q0::NamedTuple{qnames},
    m0::NamedTuple{mnames, <:Tuple{Vararg{<:TMI.MassFraction}}},
    c0::NamedTuple{names, <:Tuple{Vararg{<:Field}}},
) where {names, qnames, mnames, T <: Real}
    m = unvec(m0, mvec)
    A = watermassmatrix(m, γ)
    b = unvec(b0, bvec)
    q = unvec(q0, qvec)
    ĉ = steadyinversion(A, b, q, γ)

    Δb = b - b0
    Δc = ĉ - c0
    Δq = q - q0

    nΔb = vec(Δb)' * I * vec(Δb)
    nΔc = vec(Δc)' * I * vec(Δc)
    nΔq = vec(Δq)' * I * vec(Δq)
    nΔm = (vec(m) .- vec(m0))' * I * (vec(m) .- vec(m0))
    return nΔb + nΔc + nΔq + nΔm
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
    wetmask = trues(ngrid)
    interior = copy(wetmask)
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
    γ = Grid(axes, wetmask, interior, wrap, Δ)

    m_template = massfractions_isotropic(γ)
    mvec0 = vec(m_template)
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
    m0 = deepcopy(m_template)
    for mi in m0
        mi.fraction[wet(mi)] .*= 0.99
    end
    c0 = map(ci -> 0.95 * ci, steadyinversion(watermassmatrix(m_template, γ), b0, q0, γ))

    @testset "reverse NamedTuple steadyinversion through unvec and watermassmatrix (b, q, m)" begin
        grad_b_ad = Enzyme.make_zero(bvec0)
        grad_q_ad = Enzyme.make_zero(qvec0)
        grad_m_ad = Enzyme.make_zero(mvec0)
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Reverse),
            loss,
            Duplicated(bvec0, grad_b_ad),
            Duplicated(qvec0, grad_q_ad),
            Duplicated(mvec0, grad_m_ad),
            Const(γ),
            Const(b0),
            Const(q0),
            Const(m0),
            Const(c0),
        )

        δ = 1e-6
        grad_b_fd = finite_difference_gradient(
            bv -> loss(bv, qvec0, mvec0, γ, b0, q0, m0, c0),
            bvec0;
            δ,
        )
        grad_q_fd = finite_difference_gradient(
            qv -> loss(bvec0, qv, mvec0, γ, b0, q0, m0, c0),
            qvec0;
            δ,
        )
        grad_m_fd = finite_difference_gradient(
            mv -> loss(bvec0, qvec0, mv, γ, b0, q0, m0, c0),
            mvec0;
            δ,
        )

        @test all(isfinite.(grad_b_ad))
        @test all(isfinite.(grad_q_ad))
        @test all(isfinite.(grad_m_ad))
        @test isapprox(grad_b_ad, grad_b_fd; rtol = 1e-4, atol = 1e-8)
        @test isapprox(grad_q_ad, grad_q_fd; rtol = 1e-4, atol = 1e-8)
        @test isapprox(grad_m_ad, grad_m_fd; rtol = 1e-4, atol = 1e-8)
    end
end
