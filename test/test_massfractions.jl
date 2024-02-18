@testset "mass fractions" begin

    TMIversion = versionlist()[6] # G14 has no remote mass fractions
    #    TMIversion = "modern_90x45x33_G14"
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

    @testset "reading from matrix" begin
        m_north = TMI.massfractions_north(A, γ)

        @test maximum(m_north) < 1.0
        @test minimum(m_north) > 0.0
    end

    @testset "solve for matrix" begin

        # get observations at surface
        # set them as surface boundary condition
        y = (θ =  readfield(TMIfile, "θ", γ),
            S = readfield(TMIfile, "Sp", γ),
            δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
            P★ = preformedphosphate(TMIversion,Alu,γ),
            δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
        )

        m̃ = massfractions(y)
        Ã = watermassmatrix(m̃, γ)

        # compare m̃ and m_true (Δ̄ = 1e-14 in my case)
        for nn in keys(m̃)
            Δ = m̃[nn]
            @test maximum(m̃[nn]) ≤ 1.0
            @test minimum(m̃[nn]) ≥ 0.0
        end

        bθ = getsurfaceboundary(y.θ)
        Ãlu = lu(Ã)
        θ̃ = steadyinversion(Ã,bθ,γ)

        # compare to c.θ
        maxmisfit = 0.05
        @test Base.maximum(y.θ - θ̃) < maxmisfit
        @test Base.minimum(y.θ - θ̃) > -maxmisfit
    end
end
