using TMI
using Test

@testset "TMI.jl" begin

    TMIversion = "TMI_2010_2012_4x4x33"
    A, Alu, γ = config(TMIversion)

    ############################
    ## trackpathways
    @testset "trackpathways" begin

        @test isapprox(maximum(sum(A,dims=2)),1.0)
        @test minimum(sum(A,dims=2))> -1e-14
        
        #- define the surface patch by the bounding latitude and longitude.
        latbox = [50,60]; # 50 N -> 60 N, for example.

        # mutable due to wraparound: don't use an immutable tuple
        lonbox = [-50,0]; # 50 W -> prime meridian

        c,γ = trackpathways(TMIversion,latbox,lonbox)

        @test maximum(c[γ.wet]) ≤ 1.0
        @test minimum(c[γ.wet]) ≥ 0.0
    end
    
    #######################################
    ## example: regeneration
    @testset "regeneration" begin

        PO₄ᴿ, γ = regeneratedphosphate(TMIversion)
        @test maximum(PO₄ᴿ[γ.wet]) < 10
        @test minimum(PO₄ᴿ[γ.wet]) ≥ 0
    end
    
    ########################################
    ## example: howoceansfilled
    @testset "howoceansfilled" begin

        # are the areas and volumes consistent?
        v = cellvolume(γ)
        area = cellarea(γ)
        @test sum(0. .< v./area .< 1000.)/length(γ.I) == 1
        
        # effectively take inverse of transpose A matrix.
        dVdd = tracerinit(γ.wet); # pre-allocate c
        dVdd[γ.wet] = Alu'\v[γ.wet]

        volume = volumefilled(TMIversion)

        # volumefill positive at surface?
        @test sum(volume .< 0) == 0
        @test minimum(volume) ≤ 0.0
    end
    
    ####################################
    ## example: surfaceorigin
    @testset "surfaceorigin" begin

        # would be better to randomize location.
        xlon = 125.26; # deg E.
        xlat = -6.38;  # deg N.
        xdepth = 3000;  # meters.

        loc = (xlon,xlat,xdepth)
        δ = nearestneighbormask(loc,γ)

        @test sum(δ) > 0

        origin, γ = surfaceorigin(TMIversion,loc)
        @test isapprox(sum(filter(!isnan,origin)),1)
        @test isapprox(sum(filter(!isnan,origin)),1)
        @test minimum(filter(!isnan,origin)) ≥ 0
        @test maximum(filter(!isnan,origin)) ≤ 1
    end
    
    #########################################
    ## filterdata
    @testset "filterdata" begin

        # first guess of change to surface boundary conditions
        # ocean values are 0
        u₀ = Vector{Float64}(undef,sum(γ.wet[:,:,1]))
        fill!(u₀,zero(Float64))

        # take synthetic, noisy observations
        y, W⁻ = sample_observations(TMIversion,"θ")
        A, Alu, γ, inputfile = config(TMIversion) 

        # a first guess: observed surface boundary conditions are perfect.
        # set surface boundary condition to the observations.
        # below surface = 0 % no internal sinks or sources.
        d₀ = tracerinit(γ.wet)
        d₀[:,:,1] = y[:,:,1]

        # check gradients in misfit_gridded_data!
        fg(x) = misfit_gridded_data(x,Alu,y,d₀,W⁻,γ.wet)
        J̃₀,gJ₀ = fg(u₀)
        fg!(F,G,x) = misfit_gridded_data!(F,G,x,Alu,y,d₀,W⁻,γ.wet)

        # check with forward differences
        ϵ = 1e-5
        ii = rand(1:sum(γ.wet[:,:,1]))
        δu = copy(u₀); δu[ii] += ϵ
        ∇f_finite = (fg(δu)[1] - fg(u₀)[1])/ϵ # `[1]` to pick out cost function

        fg!(J̃₀,gJ₀,(u₀+δu)./2) # J̃₀ is not overwritten
        ∇f = gJ₀[ii]
        
        # error less than 10 percent?
        @test (∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1
        
        # sophisticated method: doesn't work due to issue #17
        #gmisfitForward = x -> ForwardDiff.gradient(misfit, x); # g = ∇f
        #gJforward = gmisfitForward(u₀)

        # gradient check with autodiff? doesn't work: issue #17
        #test = gradient(misfit,u₀)

        # filter the data with an Optim.jl method
        out = filterdata(u₀,Alu,y,d₀,W⁻,γ.wet)

        # was cost function decreased?
        @test out.minimum < J̃₀

        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J̃,gJ̃ = fg(ũ)
        @test J̃ < J̃₀

    end
    
end
