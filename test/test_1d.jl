@testset "1-dimensional domain" begin

    ngrid = (50) # number of grid cells
    xmax = 1000.0 # domain size 
    lon = collect(range(0.0,1000.0,length=ngrid[1]))
    tracer = collect(1.0.-lon./xmax)

    axes = (lon,)
    wet = trues(ngrid)
    interior = copy(wet)
    interior[begin] = false
    interior[end] = false

    wrap = (false,)
    Δ = [CartesianIndex(1,),CartesianIndex(-1,)]
    γ = Grid(axes,wet,interior,wrap,Δ)
    n = neighbors(γ)
    m0 = massfractions_isotropic(γ)
    c = Field(tracer,
        γ,
        :c,
        "linear equilibrated tracer",
        "μmol/kg")

    y = (c = c,)
    w = (c = 0.01,)

    m = massfractions(y, w)
    A = watermassmatrix(m, γ)

    ## following ex0
    # manually pick out boundary conditions
    dim =1
    b = (west = TMI.getboundarycondition(c, dim, 1, γ),
        east = TMI.getboundarycondition(c, 1, ngrid[dim], γ))

        # no real need to worry about LU decomposition in small problem.
    c̃ = steadyinversion(A,b,γ)
    @test abs(maximum(c-c̃)) < 1e-4

    # Introduce a nonconservative variable with a source
    qfield = 1.0e-2 * ones(ngrid)
    # requires negative sign which is counterintutive (needs to be fixed)
    q = TMI.Source(-qfield, γ, :q, "remineralized stuff", "μmol/kg", false)
    c_noncons = steadyinversion(A,b,γ, q=q)
    Δc = c_noncons - c
    @test iszero(sum(Δc .< 0.0))

    # ex1: trackpathways
    # no `trackpathways` function for 1D has been implemented
    # Here is a manual tracking of pathways from western boundary

     # ones() and trues(): trick to get a 0-dim array
    bwest = (west = BoundaryCondition(ones(),(),0.0,dim,1,trues()),
        east = BoundaryCondition(zeros(),(),xmax,dim,ngrid[1],trues()))

    gwest = steadyinversion(A,bwest,γ) # fraction from west boundary
    @test maximum(gwest) ≤ 1.0
    @test minimum(gwest) ≥ 0.0

    @testset "1D interpolation helpers" begin
        locvec = [lon[10]]

        δ = interpweights(locvec, γ)
        @test δ !== nothing
        @test isapprox(sum(δ), 1.0)
        @test isapprox(sum(δ .* lon), locvec[1]; atol = 1e-12)

        @test_throws ArgumentError interpweights([lon[1], lon[2]], γ)
    end
    
end 
