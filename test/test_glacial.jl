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

        # remineralized phosphate
        PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)

        # total (observed) phosphate
        PO₄total = PO₄ᴿ + PO₄pre

        ## compute total phosphate directly
        PO₄direct = steadyinversion(Alu,bPO₄,γ,q=qPO₄)
        
        ## how big is the maximum difference?
        # could replace with abs function
        @test maximum(PO₄direct - PO₄total) < 1e-2
        @test minimum(PO₄direct - PO₄total) > -1e-2

        # also compare to original field from years gone by
        @test maximum(PO₄direct - yPO₄) < 1e-2
        @test minimum(PO₄direct - yPO₄) > -1e-2

        #calculate the volume of ocean filled by surface cells in m³
        vsum = sum(10 .^ vec(volumefilled(TMIversion,Alu,γ)) .* vec(cellarea(γ)))

        #compare the volume filled to actual ocean volume
        @test isapprox(vsum, sum(cellvolume(γ)))

    end
end
