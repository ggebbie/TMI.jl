@testset "glacial" begin
    TMIversion = "LGM_90x45x33_G14"
    A, Alu, γ, TMIfile, L, B = config(TMIversion);

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
end
