@testset "mass fractions" begin
    using JuMP, COSMO
    
    TMIversion = versionlist()[6] # G14 has no remote mass fractions
    #    TMIversion = "modern_90x45x33_G14"
    A, Alu, γ, TMIfile, L, B = config(TMIversion);

    ϵ = 1e-4
    
    @testset "reading from matrix" begin
        m_north = TMI.massfractions_north(A, γ)

        @test maximum(m_north) < 1.0 + ϵ
        @test minimum(m_north) > 0.0 - ϵ
    end

    @testset "solve for matrix" begin

        for scenario in ("just determined","underdetermined")

            if scenario == "just determined"
                # get observations at surface
                # set them as surface boundary condition
                y = (θ =  readfield(TMIfile, "θ", γ),
                    S = readfield(TMIfile, "Sp", γ),
                    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
                    P★ = preformedphosphate(TMIversion,Alu,γ),
                    δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
                )

                w = (θ =  0.01,
                    S = 0.001,
                    δ¹⁸O = 0.05,
                    P★ = 0.05,
                    δ¹³C★ = 0.05
                )

            elseif scenario == "underdetermined"
            
                y = (θ =  readfield(TMIfile, "θ", γ),
                    S = readfield(TMIfile, "Sp", γ),
                )

                w = (θ =  0.01,
                    S = 0.001
                )

            end

            m̃ = massfractions(y, w)
            Ã = watermassmatrix(m̃, γ)

            # compare m̃ and m_true (Δ̄ = 1e-14 in my case)
            for nn in keys(m̃)
                Δ = m̃[nn]
                @test maximum(m̃[nn]) < 1.0 + ϵ
                @test minimum(m̃[nn]) > 0.0 - ϵ
            end

            bθ = getsurfaceboundary(y.θ)
            Ãlu = lu(Ã)
            θ̃ = steadyinversion(Ã,bθ,γ)

            # compare to c.θ
            maxmisfit = 0.2
            @test Base.maximum(y.θ - θ̃) < maxmisfit
            @test Base.minimum(y.θ - θ̃) > -maxmisfit


            # compare against the "truth" as solved by TMI
            m_true = (north = TMI.massfractions_north(A,γ),
                east   = TMI.massfractions_east(A,γ),
                south  = TMI.massfractions_south(A,γ),
                west   = TMI.massfractions_west(A,γ),
                up     = TMI.massfractions_up(A,γ),
                down   = TMI.massfractions_down(A,γ))

            #mc_test = TMI.tracer_contribution(y.θ,m_true.north) # a test
            mc_true = TMI.tracer_contribution(y.θ,m_true)
            @test maximum(mc_true) < 1e-2
            @test minimum(mc_true) > -1e-2

            mc_test = TMI.tracer_contribution(y.θ,m̃)
            @test maximum(mc_test) < 5e-2
            @test minimum(mc_test) > -5e-2
        end
    end
end
