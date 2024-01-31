@testset "mass fractions" begin
    TMIversion = "modern_90x45x33_G14"
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);
    m_north = TMI.massfractions_north(A, γ)

    @test maximum(m_north) < 1.0
    @test minimum(m_north) > 0.0
end
