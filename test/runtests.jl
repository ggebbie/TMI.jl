using Test
using TMI

function compare_controls(ubc,ubc2,testval)
    irand = rand(1:sum(ubc.wet))
    iloc = findall(ubc.wet)[irand]
    @test ubc.tracer[iloc] == ubc2.tracer[iloc]
    @test ubc.tracer[iloc] == testval
end

@testset "TMI.jl" begin

    include("test_nordic.jl")

    include("test_modern.jl")

    include("test_mitgcm.jl")
end
