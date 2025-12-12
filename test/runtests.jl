using Test
using TMI
using LinearAlgebra: lu
using SparseArrays
using Random 

percent_difference(a, b) = @. 100 * ((a - b) / b)

function centered_finite_difference(f, x; δ = 1e-3)
    orig_size = size(x)
    x_vec = vec(copy(x))
    grad_vec = similar(x_vec)
    for i in eachindex(x_vec)
        x_plus = copy(x_vec)
        x_minus = copy(x_vec)
        x_plus[i] += δ
        x_minus[i] -= δ
        grad_vec[i] = (f(x_plus) - f(x_minus)) / (2δ)
    end
    reshape(grad_vec, orig_size)
end

function compare_controls(ubc,ubc2,testval)
    irand = rand(1:sum(ubc.wet))
    iloc = findall(ubc.wet)[irand]
    @test ubc.tracer[iloc] == ubc2.tracer[iloc]
    @test ubc.tracer[iloc] == testval
end

@testset "TMI.jl" begin

    include("test_nordic.jl")

    include("test_modern.jl")
    
    include("test_glacial.jl")

    include("test_mitgcm.jl")

    include("test_massfractions.jl")

    include("test_1d.jl")

    include("test_inversion.jl")
end
