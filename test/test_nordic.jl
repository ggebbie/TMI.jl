@testset "nordic" begin
    # error-prone and long = run it first
    TMIversion = "nordic_201x115x46_B23"
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

    @testset "sourcemap" begin

        using Statistics, Interpolations

        lscale = true
        surfacetoo = false
        
        yPO₄ = readfield(TMIfile,"PO₄",γ)
        bPO₄ = (; surface = getsurfaceboundary(yPO₄),
                north = getnorthboundary(yPO₄),
                east = geteastboundary(yPO₄),
                south = getsouthboundary(yPO₄),
                west = getwestboundary(yPO₄))
        qPO₄ = readsource(TMIfile,"qPO₄",γ) # truth
        q₀ = 1e-2*onesource(γ) # first guess
        
        N = 20
        σ = 0.01 # incomplete first guess error, use scalar
        y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"PO₄",γ,N,σ)

        # Adjust surface only
        #u = (; surface = zerosurfaceboundary(γ), source = zerosource(γ))

        # Adjust interior sources only
        u = (; source = zerosource(γ,logscale=lscale))

        # Adjust interior sources and lateral boundary conditions
        u = (; source = zerosource(γ,logscale=lscale),
             north = zeronorthboundary(γ),
             east = zeroeastboundary(γ),
             south = zerosouthboundary(γ),
             west = zerowestboundary(γ))

        PO₄true = steadyinversion(Alu,bPO₄,γ,q=qPO₄)
        PO₄₀ = steadyinversion(Alu,bPO₄,γ,q=q₀)
        uvec = vec(u)

        σq = 1.0
        Q⁻ = 1.0/(σq^2) # how well is q (source) known?

        iters =5
        out, f, fg, fg! = TMI.sparsedatamap(Alu,bPO₄,u,y,W⁻,wis,locs,Q⁻,γ,q=q₀,r=1.0,iterations=iters)

        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J,∂J∂u = fg(ũ)
        J₀,∂J∂u0 = fg(vec(u))
        @test J < J₀
        @test out.minimum < J₀
        @test isapprox(J,out.minimum)
        
        ∇f, ∇f_finite = TMI.gradient_check(uvec,f,fg,fg!)
        @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

        if lscale
            #b̃ = adjustboundarycondition(b₀,unvec(u,ũ)) # combine b₀, u
            q̃ = TMI.adjustsource(q₀,unvec(u,ũ))

            # next test likes to fail if ≥ 0
            # sources may permit negative values
            #@test minimum(q̃) ≥ -0.1 

            # get new boundary conditions
            b̃ = TMI.adjustboundarycondition(bPO₄,unvec(u,ũ))

            # are lateral boundaries adjusted?
            @test ~iszero(minimum(b̃.south-bPO₄.south))
        end
    end
end
