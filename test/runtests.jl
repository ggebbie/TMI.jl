using TMI
using Test

@testset "TMI.jl" begin
    # Write your tests here.
    url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
    inputdir = "../data"

    A, Alu, γ = config(url,inputdir)
    
    c = randn(length(γ.I))
    cfld = vec2fld(c,γ.I)

    @test isapprox( sum(replace(cfld,NaN => 0.0)), sum(c))
    
    # do the functions represent a true inverse?
    @test isapprox( fld2vec(vec2fld(c,γ.I),γ.I) , c)

    # are the areas and volumes consistent?
    v = cellvolume(γ)
    area = cellarea(γ)
    @test sum(0. .< v./area .< 1000.)/length(γ.I) == 1

 
    # effectively take inverse of transpose A matrix.
    dVdd = Alu'\v

    # dVdd positive at surface?
    Isfc = surfaceIndex(γ.I)
    @test sum(dVdd[Isfc] .< 0) == 0

end
