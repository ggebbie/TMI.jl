using TMI, Test

@testset "TMI.jl" begin

    TMIversion = "TMI_2010_2012_4x4x33"
    A, Alu, γ, TMIfile = config(TMIversion)
    
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

        # choose linear or nearest neighbor interpolation
        linearinterp = true #true
        
        # would be better to randomize location.
        # xlon = 125.26; # deg E.
        # xlat = -6.38;  # deg N.
        # xdepth = 2578.2;  # meters.
        # loc = (xlon,xlat,xdepth)
        loc = wetlocation(γ)

        # choose linear or nearest neighbor interpolation

        if linearinterp
            # temporary tracer to initialize linear interp
            ctmp = tracerinit(γ.wet)
            δ = interpweights(loc,γ)
        else
            # alternate choice
            δ = nearestneighbormask(loc,γ)
        end

        @test isapprox(sum(filter(!isnan,δ)),1.0) 

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
        u₀ = zeros(Float64,sum(γ.wet[:,:,1]))

        # take synthetic, noisy observations
        y, W⁻, ctrue = sample_observations(TMIversion,"θ")

        # a first guess: observed surface boundary conditions are perfect.
        # set surface boundary condition to the observations.
        # below surface = 0 % no internal sinks or sources.
        d₀ = tracerinit(γ.wet)
        d₀[:,:,1] = y[:,:,1]

        # check gradients in misfit_gridded_data!
        fg(x) = costfunction_obs(x,Alu,d₀,y,W⁻,γ)
        f(x) = fg(x)[1]
        J̃₀,gJ₀ = fg(u₀)
        fg!(F,G,x) = costfunction_obs!(F,G,x,Alu,d₀,y,W⁻,γ)
        #fg!(J̃₀,gJ₀,u₀)
        # filter the data with an Optim.jl method

        out = filterdata(u₀,Alu,y,d₀,W⁻,fg!,γ)

        # check with forward differences
        ϵ = 1e-5
        ii = rand(1:sum(γ.wet[:,:,1]))
        δu = copy(u₀); δu[ii] += ϵ
        ∇f_finite = (f(δu) - f(u₀))/ϵ 

        fg!(J̃₀,gJ₀,(u₀+δu)./2) # J̃₀ is not overwritten
        ∇f = gJ₀[ii]
        
        # error less than 10 percent?
        @test (∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

        # was cost function decreased?
        @test out.minimum < J̃₀

        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J̃,gJ̃ = fg(ũ)
        @test J̃ < J̃₀

    end
    
    #########################################
    ## globalmap
    @testset "globalmap" begin
        # first guess of change to surface boundary conditions
        # how many randomly sampled observations?
        N = 20

        # ocean values are 0
        u₀ = zeros(Float64,sum(γ.wet[:,:,1]))

        # take synthetic, noisy observations
        y, W⁻, ctrue, locs, wis = sample_observations(TMIversion,"θ",20)

        # make a silly first guess for surface
        d₀ = tracerinit(γ.wet)
        [d₀[γ.I[ii]] = 15.0 for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]

        # check gradients in misfit_gridded_data!
        fg(x) = costfunction_obs(x,Alu,d₀,y,W⁻,wis,locs,γ)
        f(x) = fg(x)[1]
        J̃₀,gJ₀ = fg(u₀)

        ## NEXT STEP: SETUP ! FUNCTION. TEST GRADIENT. OPTIMIZE.
        #fg!(F,G,x) = costfunction_obs!(F,G,x,Alu,d₀,y,W⁻,wis,γ)
        

    end
    
end
