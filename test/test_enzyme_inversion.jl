using Enzyme, NCDatasets, Optimization, OptimizationIpopt, SparseArrays, Test, TMI
using LinearAlgebra: Diagonal, diag, issymmetric, lu
# internal functions of the extension
const TMIEnzymeOptimization = Base.get_extension(TMI, :TMIEnzymeOptimizationExt)

"""
    inversioninputs(n=5)

Small n×n×n test problem with an adjusted tracer `active` and a fixed tracer
`fixed`, each with its own interior source.
"""
function inversioninputs(n=5)
    coordinates = collect(range(-1.0, 1.0; length=n))
    interior = falses(n, n, n)
    interior[2:end-1, 2:end-1, 2:end-1] .= true
    Δ = CartesianIndex.(((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)))
    γ = Grid((180coordinates, 60coordinates, 500(coordinates .+ 1)), trues(n, n, n),
        interior, (false, false, false), collect(Δ))
    fields = (active=Field(ones(n, n, n), γ, :active, "active", "unitless"),
        fixed=Field(fill(2.0, n, n, n), γ, :fixed, "fixed", "unitless"))
    b₀ = map(getsurfaceboundary, fields)
    q₀ = (active_source=onesource(γ, :active_source, "active", "unitless"),
        fixed_source=2onesource(γ, :fixed_source, "fixed", "unitless"))
    r = (active=(active_source=1.0,), fixed=(fixed_source=1.0,))
    m₀ = massfractions_isotropic(γ)
    foreach(m -> replace!(m.fraction, NaN => 0.0), m₀)
    σm = unvec(m₀, collect(range(0.2, 0.8; length=length(vec(m₀)))))
    c = steadyinversion(lu(watermassmatrix(m₀, γ)), b₀, q₀, r, γ)
    y = Observations(map(cₖ -> 0.95cₖ, c), 0.5, γ)
    controls = (b=(active=(σ=0.5,),),
        q=(active_source=(σ=1.0, lb=0.25, ub=4.0, logscale=true),),
        m=(σ=σm, lb=0.0, ub=1.0))
    inversion = Inversion(y; b₀, q₀, r, m₀, controls)
    locs = [(γ.lon[i], γ.lat[i], γ.depth[i]) for i in (2, 2, 4)]
    return (; γ, b₀, q₀, r, m₀, σm, c, y, controls, inversion, locs)
end

@testset "Enzyme inversion" begin
    (; γ, b₀, q₀, r, m₀, σm, c, y, controls, inversion, locs) = inversioninputs()
    ninterior = sum(γ.interior)

    @testset "Observations" begin
        @test vec(y) == vcat(y.tracers.active.tracer[γ.interior], y.tracers.fixed.tracer[γ.interior])
        @test diag(y.Wⁱ) ≈ fill(0.5^-2, 2ninterior)
        point_values = observe(c.fixed, locs, γ)
        mixed = Observations((active=y.tracers.active, fixed=point_values), (active=0.5, fixed=fill(0.5, 3)), γ;
            locs=(fixed=locs,))
        @test length(vec(mixed)) == ninterior + 3
        σ = Field(fill(2.0, size(γ.wet)), γ, :σ, "uncertainty", "unitless")
        L = Field(fill(1_000.0, size(γ.wet)), γ, :L, "scale", "km")
        correlated = Observations((active=c.active,), σ, γ; L)
        @test correlated.Wⁱ isa SparseMatrixCSC
        @test issymmetric(correlated.Wⁱ)
        given = Observations((active=point_values,), 0.5, γ; locs=locs, Wⁱ=Diagonal(fill(4.0, 3)))
        @test given.Wⁱ ≈ Diagonal(fill(4.0, 3))
        @test_throws ArgumentError Observations((active=point_values,), 0.5, γ; locs=locs, L=1_000.0)
    end

    @testset "Inversion" begin
        @test inversion.ranges.b == 1:length(vec(b₀.active))
        @test length(inversion.ranges.q) == ninterior
        @test length(inversion.ranges.m) == length(vec(m₀))
        uvec = zero(inversion.σ)
        b, q, m = TMIEnzymeOptimization.adjustfirstguess!(deepcopy(inversion.u), inversion, uvec)
        @test vec(b) == vec(b₀)
        @test vec(q.active_source) ≈ vec(q₀.active_source)
        @test vec(m) == vec(m₀)
        uvec[inversion.ranges.q] .= log(3.0)
        _, q, _ = TMIEnzymeOptimization.adjustfirstguess!(deepcopy(inversion.u), inversion, uvec)
        @test vec(q.active_source) ≈ fill(3.0, ninterior)
        @test vec(q.fixed_source) == vec(q₀.fixed_source)
        @test inversion.lb[inversion.ranges.q] ≈ fill(log(0.25), ninterior)
        @test inversion.ub[inversion.ranges.q] ≈ fill(log(4.0), ninterior)
        @test inversion.lb[inversion.ranges.m] ≈ -vec(m₀)
        @test inversion.σ[inversion.ranges.m] ≈ vec(σm)
        linear = Inversion(y; b₀, q₀, r, m₀, controls=(q=(active_source=(σ=1.0, lb=-4.0, ub=4.0),),))
        @test linear.lb[linear.ranges.q] ≈ fill(-5.0, ninterior)
        @test isempty(linear.ranges.b)
        @test isempty(linear.ranges.m)
    end

    @testset "Cost" begin
        uvec = 0.01 .* inversion.σ
        cache = TMIEnzymeOptimization.EnzymeCostGradientCache(inversion)
        terms = TMIEnzymeOptimization.costterms(
            TMIEnzymeOptimization.steadyinversion(uvec, inversion, cache.u, cache.F), uvec, inversion)
        @test terms.Jcontrol ≈ sum(abs2, uvec ./ inversion.σ)
        # cost alone, as for Ipopt's trial points, against the primal of the Enzyme pass
        J = TMIEnzymeOptimization.cost!(cache, uvec)
        @test J ≈ terms.Jdata + terms.Jcontrol
        TMIEnzymeOptimization.costgradient!(cache, uvec)
        @test cache.J ≈ J
        smooth = Inversion(y; b₀, q₀, r, m₀, controls=(b=(active=(σ=0.5, L=1_000.0),),))
        @test smooth.Q⁻ ≈ smooth.Q⁻'
        @test all(diag(smooth.Q⁻) .≥ smooth.σ .^ -2)
    end

    # gradient_check errors when Enzyme and finite differences disagree
    checkgradient(inversion) = TMI.gradient_check(inversion; n=1) isa NamedTuple
    @testset "Gradient" begin
        @test checkgradient(inversion)
        # the gradient where F already holds A, so lu! is skipped, against one that refactors F
        uvec = 0.01 .* inversion.σ
        skipped, refactored = (TMIEnzymeOptimization.EnzymeCostGradientCache(inversion) for _ in 1:2)
        TMIEnzymeOptimization.cost!(skipped, uvec)
        TMIEnzymeOptimization.costgradient!(skipped, uvec)
        TMIEnzymeOptimization.costgradient!(refactored, uvec)
        @test skipped.guvec == refactored.guvec
        mixed = Observations((active=y.tracers.active, fixed=observe(c.fixed, locs, γ)), 0.5, γ;
            locs=(fixed=locs,))
        @test checkgradient(Inversion(mixed; b₀, q₀, r, m₀, controls))
        @test checkgradient(Inversion(y; b₀, q₀, r, m₀, controls=(b=controls.b, q=controls.q)))
        @test checkgradient(Inversion(y; b₀, q₀, r, m₀, controls=(m=controls.m,)))
        # two adjusted boundary conditions and two adjusted sources with fixed m;
        # Enzyme failed to compile this on Julia 1.13 before `unvec!` stopped looping over NamedTuples
        @test checkgradient(Inversion(y; b₀, q₀, r, m₀, controls=(b=(active=(σ=0.5,), fixed=(σ=0.5,)),
            q=(active_source=controls.q.active_source, fixed_source=(σ=1.0,)))))
    end

    @testset "Mass conservation" begin
        conservation = TMIEnzymeOptimization.massconservation(inversion)
        z = zero(inversion.σ)
        values = zeros(ninterior)
        conservation.cons(values, z, nothing)
        @test values ≈ ones(ninterior)
        i = first(inversion.ranges.m)
        δ = zero(z)
        δ[i] = 1.0e-6
        plus, minus = similar(values), similar(values)
        conservation.cons(plus, z + δ, nothing)
        conservation.cons(minus, z - δ, nothing)
        @test (plus - minus) / 2.0e-6 ≈ Vector(conservation.cons_jac_prototype[:, i])
        @test isnothing(TMIEnzymeOptimization.massconservation(
            Inversion(y; b₀, q₀, r, m₀, controls=(b=controls.b,))).cons)
    end

    @testset "runinversion" begin
        mixed = Observations((active=y.tracers.active, fixed=observe(c.fixed, locs, γ)), 0.5, γ;
            locs=(fixed=locs,))
        mixed_inversion = Inversion(mixed; b₀, q₀, r, m₀, controls)
        mktempdir() do directory
            checkpoint_directory = joinpath(directory, "checkpoints")
            data_directory = joinpath(directory, "data")
            runinversion(mixed_inversion, IpoptOptimizer(hessian_approximation="limited-memory");
                name="test", iterations=2, checkpoint_interval=2,
                checkpoint_directory, data_directory, number_of_gradient_checks=0)
            NCDataset(joinpath(checkpoint_directory, "test_observations.nc")) do ds
                @test ds["fixed_sample"][:] ≈ observe(c.fixed, locs, γ)
                @test ds["fixed"].attrib["support"] == "sparse"
                @test ds["active"].attrib["support"] == "gridded"
            end
            NCDataset(joinpath(checkpoint_directory, "test_iteration_000000.nc")) do ds
                @test ds.attrib["iteration"] == 0
                cache = TMIEnzymeOptimization.EnzymeCostGradientCache(mixed_inversion)
                @test ds.attrib["J"] ≈ TMIEnzymeOptimization.cost!(cache, zero(mixed_inversion.σ))
                @test ds.attrib["J"] ≈ ds.attrib["Jdata"] + ds.attrib["Jcontrol"]
                @test haskey(ds, "active_source")
            end
            @test isfile(joinpath(checkpoint_directory, "test_iteration_000002.nc"))
            @test isfile(joinpath(checkpoint_directory, "ipopt_output.txt"))
            @test isfile(joinpath(data_directory, "test_final.nc"))
        end
    end
end
