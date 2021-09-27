using TMI
using Test

@testset "TMI.jl" begin
    # Write your tests here.
    
    c = randn(length(γ.I))
    cfld = vec2fld(c,γ.I)

    @test isapprox( sum(replace(cfld,NaN => 0.0)), sum(c))
    
    # do the functions represent a true inverse?
    @test isapprox( fld2vec(vec2fld(c,γ.I),γ.I) , c)
end
