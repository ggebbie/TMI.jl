using TMI
using Test

@testset "TMI.jl" begin

    ############################
    ## preamble
    url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
    inputdir = "../data"
    
    A, Alu, c, ΔPO₄, γ = config(url,inputdir)

    PO₄ᴿ = tracerinit(γ.wet); # pre-allocate 
    PO₄ᴿ[γ.wet] = Alu\ΔPO₄[γ.wet]

    @test minimum(ΔPO₄[γ.wet]) ≥ 0
    @test isapprox(maximum(sum(A,dims=2)),1.0)
    @test minimum(sum(A,dims=2))> -1e-14

    #################################
    ## example: trackpathways
    
    #- define the surface patch by the bounding latitude and longitude.
    latbox = [50,60]; # 50 N -> 60 N, for example.

    # mutable due to wraparound: don't use an immutable tuple
    lonbox = [-50,0]; # 50 W -> prime meridian

    d = surfacepatch(lonbox,latbox,γ)

    # do matrix inversion to get quantity of dyed water throughout ocean:
    c = tracerinit(γ.wet); # pre-allocate c

    # make methods that make the "wet" index unnecessary
    c[γ.wet] = Alu\d[γ.wet] # presumably equivalent but faster than `c = A\d`

    @test maximum(c[γ.wet]) ≤ 1.0
    @test minimum(c[γ.wet]) ≥ 0.0

    #######################################
    ## example: regeneration

    PO₄ᴿ = tracerinit(γ.wet); # pre-allocate 
    PO₄ᴿ[γ.wet] = -(Alu\ΔPO₄[γ.wet])
    @test maximum(PO₄ᴿ[γ.wet]) < 10
    @test minimum(PO₄ᴿ[γ.wet]) ≥ 0

    ########################################
    ## example: howoceansfilled
    
    # are the areas and volumes consistent?
    v = cellvolume(γ)
    area = cellarea(γ)
    @test sum(0. .< v./area .< 1000.)/length(γ.I) == 1
 
    # effectively take inverse of transpose A matrix.
    dVdd = tracerinit(γ.wet); # pre-allocate c
    dVdd[γ.wet] = Alu'\v[γ.wet]

    # scale the sensitivity value by surface area so that converging meridians are taken into account.
    I = γ.I
    volumefill = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    fill!(volumefill,0.0)

    [volumefill[I[ii][1],I[ii][2]] = dVdd[I[ii]] / area[I[ii][1],I[ii][2]] for ii ∈ eachindex(I) if I[ii][3] == 1]

    # volumefill positive at surface?
    @test sum(volumefill .< 0) == 0
    @test minimum(volumefill) ≤ 0.0

    ####################################
    ## example: surfaceorigin

    # would be better to randomize location.
    xlon = 125.26; # deg E.
    xlat = -6.38;  # deg N.
    xdepth = 3000;  # meters.
    loc = (xlon,xlat,xdepth)

    δ = nearestneighbormask(loc,γ)
    @test isapprox(sum(replace(δ,NaN=>0.0)),1)

    dVdδ = tracerinit(γ.wet); # pre-allocate c
    dVdδ[γ.wet] = Alu'\δ[γ.wet]

    # surfaceorigin only exists at sea surface
    surfaceorigin = view(dVdδ,:,:,1)

    @test isapprox(sum(filter(!isnan,surfaceorigin)),1)
    @test minimum(filter(!isnan,surfaceorigin)) ≥ 0
    @test maximum(filter(!isnan,surfaceorigin)) ≤ 1
    
end
