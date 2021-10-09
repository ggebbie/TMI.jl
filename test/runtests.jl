using TMI
using Test

@testset "TMI.jl" begin

    ## example 1
    url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
    inputdir = "../data"
    
    Azyx, Axyz, Alu, c, γ = config(url,inputdir)

    @test isapprox(maximum(sum(Axyz,dims=2)),1.0)
    @test minimum(sum(Axyz,dims=2))> -1e-14

    #- define the surface patch by the bounding latitude and longitude.
    latbox = [50,60]; # 50 N -> 60 N, for example.

    # mutable due to wraparound: don't use an immutable tuple
    lonbox = [-50,0]; # 50 W -> prime meridian

    d = surfacePatch(lonbox,latbox,γ)

    # do matrix inversion to get quantity of dyed water throughout ocean:
    c = tracerFieldInit(γ.wet); # pre-allocate c

    # make methods that make the "wet" index unnecessary
    c[γ.wet] = Alu\d[γ.wet] # presumably equivalent but faster than `c = A\d`

    @test (maximum(c[γ.wet]) ≤ 1.0)
    @test (minimum(c[γ.wet]) ≥ 0.0)
    
    ## example 2
    # are the areas and volumes consistent?
    v = cellVolume(γ)
    area = cellArea(γ)
    @test sum(0. .< v./area .< 1000.)/length(γ.I) == 1
 
    # effectively take inverse of transpose A matrix.
    dVdd = tracerFieldInit(γ.wet); # pre-allocate c
    dVdd[γ.wet] = Alu'\v[γ.wet]

    # scale the sensitivity value by surface area so that converging meridians are taken into account.
    I = γ.I
    volumefill = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    fill!(volumefill,0.0)

    [volumefill[I[ii][1],I[ii][2]] = dVdd[I[ii]] / area[I[ii][1],I[ii][2]] for ii ∈ eachindex(I) if I[ii][3] == 1]

    # volumefill positive at surface?
    @test sum(volumefill .< 0) == 0
    @test minimum(volumefill) ≤ 0.0

end
