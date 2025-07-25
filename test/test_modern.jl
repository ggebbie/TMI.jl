@testset "modern" begin
    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_180x90x33_GH11_GH12"
    #TMIversion = "modern_90x45x33_unpub12"
    
    A, Alu, γ, TMIfile, L, B = config(TMIversion);

    @testset "trim domain" begin

        b_south = ones(2,45,γ,:bc_south,"South Boundary","nondim")
        b_up = ones(3,25,γ,:bc_upper,"Upper Boundary","nondim")
        b_lo = ones(3,29,γ,:bc_lower,"Lower Boundary","nondim")
        b_surface = ones(3,1,γ,:bc_surface,"Surface","nondim")

        # find all locations where boundary set to one
        b_ones = (surface = b_surface,
          south = b_south,
          upper = b_up,
          lower = b_lo)

        # update grid to be consistent with boundary conditions
        γb = Grid(b_ones, γ)
        @test sum(γ.interior) > sum(γb.interior) # interior points replaced with boundaries
        
        # update water-mass matrix to be consistent with boundary points and grid 
        Abc = watermassmatrix(A, γb) # add test here
    end

    @testset "regions" begin
        # test that all global values sum to 1.
        Nr = length(TMI.regionlist())
        sumb = zeros(Int,Nr)
        for i in 1:Nr
            # read v1 of regions from NetCDF file: used Floating point numbers for mask
            sumb[i] = sum(TMI.surfaceregion(TMIversion,TMI.regionlist()[i]).tracer)
        end

        @test sum(sumb[1]) == sum(sumb[2:8])
        @test sum(sumb[1]) == sum(sumb[vcat(5,7:Nr)])

    end

    @testset "steadyinversion" begin

        # get observations at surface
        # set them as surface boundary condition
        yPO₄ = readfield(TMIfile,"PO₄",γ)

        # a first guess: observed surface boundary conditions are perfect.
        # set surface boundary condition to the observations.

        # choice: BoundaryCondition or NamedTuple(BoundaryCondition)
        #bPO₄ = getsurfaceboundary(yPO₄)
        bPO₄ = getsurfaceboundary(yPO₄)

        ## preformed phosphate
        PO₄pre = steadyinversion(Alu,bPO₄,γ)

        ## read phosphate source
        qPO₄ = readsource(TMIfile,"qPO₄",γ)

        # zero boundary condition, choose one line of next two
        #b₀ = zerosurfaceboundary(γ)
        b₀ = zerosurfaceboundary(γ)

        # side check: onesurfaceboundary retains units
        b_argon39 = onesurfaceboundary(γ,:Ar39,"argon-39","% modern")
        @test b_argon39.name == :Ar39
        
        # remineralized phosphate
        PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)

        # total (observed) phosphate
        PO₄total = PO₄ᴿ + PO₄pre

        ## compute total phosphate directly
        PO₄direct = steadyinversion(Alu,bPO₄,γ,q=qPO₄)
        
        ## how big is the maximum difference?
        # could replace with abs function
        @test maximum(PO₄direct - PO₄total) < 0.1
        @test minimum(PO₄direct - PO₄total) > -0.1

        ## oxygen distribution, just be sure it runs
        yO₂ = readfield(TMIfile,"O₂",γ)
        bO₂ = getsurfaceboundary(yO₂)
        O₂ = steadyinversion(Alu,bO₂,γ,q=qPO₄,r=-170.0)

        ## radio-decay
        Aradio = watermassmatrix(TMIfile, γ, 5730) # like radiocarbon
        Oradio = steadyinversion(Aradio,bO₂,γ)

        Aradio2 = watermassmatrix(TMIfile, γ, 269) # like radiocarbon
        Oradio2 = steadyinversion(Aradio2,bO₂,γ)

        @test minimum(Oradio) > minimum(Oradio2)
    end

    ############################
    ## trackpathways
    @testset "trackpathways" begin

        @test isapprox(maximum(sum(A,dims=2)),1.0)
        @test minimum(sum(A,dims=2))> -1e-14
        
        latbox = [50,60]; # 50 N -> 60 N, for example.
        lonbox = [-50, 0]; # 50 W -> prime meridian

        #b = surfacepatch(lonbox,latbox,γ)
        c = trackpathways(Alu,latbox,lonbox,γ)

        @test maximum(c) ≤ 1.0
        @test minimum(c) ≥ 0.0
    end

    @testset "watermassdistribution" begin
        list = TMI.regionlist()
        region = list[2]

        #b = TMI.surfaceregion(TMIversion,region,γ)
        #g = steadyinversion(Alu,b,γ)
        g = watermassdistribution(TMIversion,Alu,region,γ)

        @test maximum(g) ≤ 1.0
        @test minimum(g) ≥ 0.0

    end
    
    #######################################
    ## example: regeneration
    @testset "regeneration" begin

        PO₄ᴿ = regeneratedphosphate(TMIversion,Alu,γ)
        @test maximum(PO₄ᴿ) < 5
        @test minimum(PO₄ᴿ) ≥ 0
    end

    #######################################
    ## example: mean age
    @testset "meanage" begin

        a = meanage(TMIversion,Alu,γ)
        @test maximum(a) < 3000
        @test minimum(a) ≥ 0
    end

    ########################################
    ## example: howoceansfilled
    @testset "howoceansfilled" begin

        # are the areas and volumes consistent?
        v = cellvolume(γ)
        area = cellarea(γ)
        #@test sum(0. .< v./area .< 1000.)/length(γ.I) == 1
        volume = volumefilled(TMIversion,Alu,γ)
        # volumefill no smaller than smallest box?
        @test exp10(minimum(volume)) ≥ 4.9
    end

    ####################################
    ## example: effective endmembers
    @testset "effective endmembers" begin
        region = TMI.regionlist()[2]
        θ = readfield(TMIfile,"θ",γ)

        @test isapprox(TMI.effective_endmember(TMIversion,Alu,θ,"ANT",γ), -1.31, atol=0.2)
        @test isapprox(TMI.effective_endmember(TMIversion,Alu,θ,"NATL",γ), 2.64, atol=0.2)
    end
    
    ####################################
    ## example: surfaceorigin
    @testset "surfaceorigin" begin

        # choose linear or nearest neighbor interpolation
        linearinterp = true #true
        # randomized location.
        loc = wetlocation(γ)

        # choose linear or nearest neighbor interpolation

        if linearinterp
            # temporary tracer to initialize linear interp
            ctmp = zeros(γ.wet)
            δ = interpweights(loc,γ)
        else
            # alternate choice
            δ = nearestneighbormask(loc,γ)
        end
        @test isapprox(sum(filter(!isnan,δ)),1.0) 

        origin = surfaceorigin(loc, Alu, γ)

        @test 0.99 < sum(exp10.(origin.tracer[origin.wet])) < 1.01
        @test maximum(origin) ≤ 0 # log10(1) = 0
    end
    
    #########################################
    ## formerly filterdata
    @testset "steadyclimatology" begin

        for ii = 1:2
            y, W⁻, ctrue = synthetic_observations(TMIversion,"θ",γ)
            if ii == 1
                println("NamedTuple type")
                u = (;surface = zerosurfaceboundary(γ))
                b = (;surface = getsurfaceboundary(y))
            else
                println("BoundaryCondition type")
                u = zerosurfaceboundary(γ)
                b = getsurfaceboundary(y)
            end

            out, f, fg, fg! = steadyclimatology(Alu,b,u,y,W⁻,γ)

            # reconstruct by hand to double-check.
            ũ = out.minimizer
            J,∂J∂u = fg(ũ)
            J₀,∂J∂u0 = fg(vec(u))
            @test J < J₀
            @test out.minimum < J₀
            @test isapprox(J,out.minimum)
                
            ∇f, ∇f_finite = TMI.gradient_check(vec(u),f,fg,fg!)
            @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

        end
    end
    
    #########################################
    ## globalmap
    @testset "sparsedatamap" begin

        using Interpolations, LinearAlgebra
        using Statistics
        
        N = 20
        y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N)

        for ii = 1:2
            if ii == 1
                println("NamedTuple type")
                u = (;surface = zerosurfaceboundary(γ))
                b = (;surface = mean(y) * onesurfaceboundary(γ))
            else
                println("BoundaryCondition type")
                u = zerosurfaceboundary(γ)
                b = mean(y) * onesurfaceboundary(γ)
            end
            uvec = vec(u)
            σb = 5.0

            for jj = 1:2
                if jj == 1
                    global Q⁻ = 1.0/(σb^2) # a scalar
                    # spatially uniform first-guess expected size
                elseif jj==2
                    lengthscale = 1000.0
                    Dg = gaussiandistancematrix(γ,σb,lengthscale)
                    Q⁻ = inv(cholesky(Dg))
                end

                iters =5
                out, f, fg, fg! = TMI.sparsedatamap(Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,iterations=iters)

                # reconstruct by hand to double-check.
                ũ = out.minimizer
                J,∂J∂u = fg(ũ)
                J₀,∂J∂u0 = fg(vec(u))
                @test J < J₀
                @test out.minimum < J₀
                @test isapprox(J,out.minimum)
                
                ∇f, ∇f_finite = TMI.gradient_check(uvec,f,fg,fg!)
                @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1
            end
        end
    end

    @testset "sourcemap" begin

        using Interpolations
        using Statistics

        for lscale in (false,true)
            for surfacetoo in (false,true)
                
                yPO₄ = readfield(TMIfile,"PO₄",γ)
                bPO₄ = getsurfaceboundary(yPO₄)
                qPO₄ = readsource(TMIfile,"qPO₄",γ)
                q₀ = 1e-2*onesource(γ)
                
                N = 20
                y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"PO₄",γ,N)

                #u = (; source = zerosource(γ))
                if surfacetoo
                    u = (; surface = zerosurfaceboundary(γ), source = zerosource(γ))
                else
                    u = (; source = zerosource(γ,logscale=lscale))
                end
                
                b = (; surface = bPO₄) # surface boundary condition

                PO₄true = steadyinversion(Alu,b,γ,q=qPO₄)
                PO₄₀ = steadyinversion(Alu,b,γ,q=q₀)
                uvec = vec(u)

                σq = 1.0
                Q⁻ = 1.0/(σq^2) # how well is q (source) known?

                iters =5
                out, f, fg, fg! = TMI.sparsedatamap(Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=q₀,r=1.0,iterations=iters)

                                # reconstruct by hand to double-check.
                ũ = out.minimizer
                J,∂J∂u = fg(ũ)
                J₀,∂J∂u0 = fg(vec(u))
                @test J < J₀
                @test out.minimum < J₀
                @test isapprox(J,out.minimum)
                
                ∇f, ∇f_finite = TMI.gradient_check(uvec,f,fg,fg!)
                @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

                if lscale
                    #b̃ = adjustboundarycondition(b₀,unvec(u,ũ)) # combine b₀, u
                    q̃ = TMI.adjustsource(q₀,unvec(u,ũ))
                    #@test minimum(q̃) ≥ 0 # likes to fail from time to time
                end
            end
        end
    end

    @testset "unvec" begin
        u = (;surface = onesurfaceboundary(γ),
             west = 2 * onewestboundary(γ),
             east = 3 *oneeastboundary(γ))

        vecu = vec(u)
        
        u2 = (;surface = zerosurfaceboundary(γ),
             west = 2 * zerowestboundary(γ),
             east = 3 * zeroeastboundary(γ))

        u3 = unvec(u2,vecu)
        compare_controls(u.surface,u3.surface,1.0)
        compare_controls(u.east,u3.east,3.0)
        compare_controls(u.west,u3.west,2.0)
        
        unvec!(u2,vecu)
        compare_controls(u.surface,u2.surface,1.0)
        compare_controls(u.east,u2.east,3.0)
        compare_controls(u.west,u2.west,2.0)

    end

    @testset "mixed layer" begin
        mixedlayer = mixedlayermask(A,γ)
        @test sum(mixedlayer) == 13744
        @test iszero(sum(mixedlayer[:,:,1]))

        τ = 0.1
        Lmix = mixedlayermatrix(A, γ, τ)
        @test Lmix isa SparseMatrixCSC
        @test maximum(Lmix) < 2/τ

        Ldir = dirichletmatrix(γ, τ)
        @test maximum(Ldir) < 2 / τ
    end
end
