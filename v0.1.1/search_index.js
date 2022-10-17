var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TMI","category":"page"},{"location":"#TMI","page":"Home","title":"TMI","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TMI.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TMI]","category":"page"},{"location":"#TMI.cartesianindexXYZ-Tuple{Any}","page":"Home","title":"TMI.cartesianindexXYZ","text":"function cartesianindexXYZ(wet)\nRead and assemble the grid coordinates\naccording to a 3D tracer in x,y,z order\n\nArguments\n\nwet: BitArray logical mask for wet points\n\nOutput\n\nI: 3D Cartesian indices\n\n\n\n\n\n","category":"method"},{"location":"#TMI.cartesianindexZYX-Tuple{Any}","page":"Home","title":"TMI.cartesianindexZYX","text":"function cartesianindexZYX(file)\nRead and assemble the grid coordinates\naccording to the legacy MATLAB code (z,y,x order).\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\ngrid: TMI grid coordinates\n\n\n\n\n\n","category":"method"},{"location":"#TMI.config-Tuple{Any}","page":"Home","title":"TMI.config","text":"function config(url,inputdir)\n\nArguments\n\nTMIversion: TMI version for water-mass/circulation model\n\nOutput\n\nA: TMI steady-state water-mass matrix\nAlu: LU decomposition of A\nγ: TMI grid properties\nTMIfile: TMI file name\n\n\n\n\n\n","category":"method"},{"location":"#TMI.control2state!-Union{Tuple{T}, Tuple{Array{T, 3}, Vector{T}, Any}} where T<:Real","page":"Home","title":"TMI.control2state!","text":"function control2state!(c,u,γ)\nAdd surface control vector to existing 3D field\n\nArguments\n\nc:: state field, 3d tracer field with NaN on dry points, modified by function\nu:: surface control vector\nwet::BitArray mask of ocean points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.control2state!-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, Any}} where T<:Real","page":"Home","title":"TMI.control2state!","text":"function control2state!(c,u,γ)\nAdd surface control vector to existing 3D field\n\nArguments\n\nc:: state field, 3d tracer field with NaN on dry points, modified by function\nu:: surface control vector\nwet::BitArray mask of ocean points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.control2state-Union{Tuple{T}, Tuple{Matrix{T}, Any}} where T<:Real","page":"Home","title":"TMI.control2state","text":"function control2state(tracer2D,γ)\nturn 2D surface field into 3D field with zeroes below surface\n\nArguments\n\ntracer2D:: 2D surface tracer field\nwet::BitArray mask of ocean points\n\nOutput\n\ntracer3D:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.control2state-Union{Tuple{T}, Tuple{Vector{T}, Any}} where T<:Real","page":"Home","title":"TMI.control2state","text":"function control2state(u,γ)\nturn surface control vector into 3D field with zeroes below surface\n\nArguments\n\nu:: surface control vector\nwet::BitArray mask of ocean points\n\nOutput\n\ntracer3D:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction!-Union{Tuple{T}, Tuple{Any, Any, Vector{T}, Any, Array{T, 3}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, Any, TMI.grid}} where T<:Real","page":"Home","title":"TMI.costfunction!","text":"function costfunction!(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,Q⁻,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\nIssue: couldn't figure out how to nest with costfunction_obs!\n\nArguments\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ndfld: model constraints\ny: pointwise observations\nWⁱ: inverse of W weighting matrix for observations\nwis: weights for interpolation (data sampling, E)\nlocs: data locations (lon,lat,depth)\nQ⁻: weights for control vector\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction-Union{Tuple{T}, Tuple{Vector{T}, Any, Array{T, 3}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, Any, TMI.grid}} where T<:Real","page":"Home","title":"TMI.costfunction","text":"function costfunction(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,Q⁻,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\nIssue: couldn't figure out how to nest with costfunction_obs!\n\nArguments\n\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ndfld: model constraints\ny: pointwise observations\nWⁱ: inverse of W weighting matrix for observations\nwis: weights for interpolation (data sampling, E)\nlocs: data locations (lon,lat,depth)\nQ⁻: weights for control vector\nγ: grid\n\nOutput\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs!-Union{Tuple{T}, Tuple{Any, Any, Vector{T}, Any, Array{T, 3}, Array{T, 3}, LinearAlgebra.Diagonal{T, Vector{T}}, TMI.grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs!","text":"function costfunction_obs(u,Alu,y,d,Wⁱ,wet)\nsquared model-data misfit for gridded data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\nd: model constraints\ny: observations on grid\nWⁱ: inverse of W weighting matrix for observations\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs!-Union{Tuple{T}, Tuple{Any, Any, Vector{T}, Any, Array{T, 3}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, TMI.grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs!","text":"function costfunction_obs!(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,locs,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ndfld: model constraints\ny: pointwise observations\nWⁱ: inverse of W weighting matrix for observations\nwis: weights for interpolation (data sampling, E)\nlocs: data locations (lon,lat,depth)\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs-Union{Tuple{T}, Tuple{Vector{T}, Any, Array{T, 3}, Array{T, 3}, LinearAlgebra.Diagonal{T, Vector{T}}, TMI.grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs","text":"function costfunction_obs(u,Alu,dfld,yfld,Wⁱ,γ)\nsquared model-data misfit for gridded data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ny: observations on grid\nd: model constraints\nWⁱ: inverse of W weighting matrix for observations\nwet: BitArray ocean mask\n\nOutput\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs-Union{Tuple{T}, Tuple{Vector{T}, Any, Array{T, 3}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, TMI.grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs","text":"function costfunction_obs(u,Alu,dfld,yfld,Wⁱ,wis,locs,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ny: pointwise observations\nd: model constraints\nWⁱ: inverse of W weighting matrix for observations\nwis: weights for interpolation \nlocs: data locations\nγ: grid\n\nOutput\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\n\n\n\n\n\n","category":"method"},{"location":"#TMI.depthindex-Tuple{Any}","page":"Home","title":"TMI.depthindex","text":"function depthindex(I) \nGet the k-index (depth level) from the Cartesian index\n\n\n\n\n\n","category":"method"},{"location":"#TMI.download-Tuple{Any, Any}","page":"Home","title":"TMI.download","text":"function download(url,inputdir)\nRead and assemble all TMI inputs.\n\nArguments\n\nurl: Google Drive location of TMI input\ninputdir: input directory location to store file\n\nOutput\n\nnone\n\n\n\n\n\n","category":"method"},{"location":"#TMI.dyeplot-NTuple{4, Any}","page":"Home","title":"TMI.dyeplot","text":"function dyeplot\nPlot of dye in ocean\n\nArguments\n\nlat: latitude arrays\ndepth: depth array\nvals: lat x depth value array\nlims: contour levels\n\n\n\n\n\n","category":"method"},{"location":"#TMI.filterdata-NTuple{7, Any}","page":"Home","title":"TMI.filterdata","text":"function filterdata(u₀,Alu,y,d₀,W⁻,γ)      Find the distribution of a tracer given:      (a) the pathways described by A or its LU decomposition Alu,      (b) first-guess boundary conditions and interior sources given by d₀,      (c) perturbations to the surface boundary condition u₀     that best fits observations, y,     according to the cost function,     J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)     subject to Aỹ = d₀ + Γ u₀.                      W⁻ is a (sparse) weighting matrix.     See Supplementary Section 2, Gebbie & Huybers 2011.\n\nArguments\n\nu₀:\nAlu:\nd₀: first guess of boundary conditions and interior sources\ny: observations on 3D grid\nW⁻: weighting matrix best chosen as inverse error covariance matrix\nfg!: compute cost function and gradient in place\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.fld2vec-Tuple{Array{Float64, 3}, Vector{CartesianIndex{3}}}","page":"Home","title":"TMI.fld2vec","text":"function fld2vec\nTransfer 3D field with accounting for ocean bathymetry to a vector without land points\n\nArguments\n\nfield: field in 3d form including land points (NaN)\nI: cartesian indices of ocean points\n\nOutput\n\nvector: field in vector form (no land points)\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gdriveurl-Tuple{Any}","page":"Home","title":"TMI.gdriveurl","text":"function gdriveurl(TMIversion)\nplaceholder function to give location (URL) of Google Drive input\nin the future, consider a struct or Dict that describes all TMI versions.\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nurl: location (URL) for download\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gridprops-Tuple{Any}","page":"Home","title":"TMI.gridprops","text":"function gridprops(file)\nRead and assemble the grid properties.\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\ngrid: TMI grid coordinates\n\n\n\n\n\n","category":"method"},{"location":"#TMI.horizontaldistance-Tuple{Any, TMI.grid}","page":"Home","title":"TMI.horizontaldistance","text":"function horizontaldistance(loc,γ)\nreturn the Cartesian index and linear index \nof the nearest N neighbors\n\nArguments\n\nloc: 3-tuple of lon,lat,depth location\nγ: TMI.grid\n\nOutput\n\nhordist: horizontal distance to nearest tracer grid points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.interpindex-Tuple{Any, Any}","page":"Home","title":"TMI.interpindex","text":"function interpindex(loc,γ)     Weights for linear interpolation.     The derivative of linear interpolation is needed in sensitivity studies.     ReverseDiff.jl could find this quantity automatically.     Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.\n\nArguments\n\nc: a temporary tracer field, would be nice to make it unnecessary\nloc: (lon,lat,depth) tuple of a location of interest\nγ: TMI grid\n\nOutput\n\nδ: weights on a 3D tracer field grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.interpweights-Tuple{Any, Any}","page":"Home","title":"TMI.interpweights","text":"function interpweights(loc,γ)     Weights for linear interpolation.     The derivative of linear interpolation is needed in sensitivity studies.     ReverseDiff.jl could find this quantity automatically.     Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.\n\nArguments\n\nloc: (lon,lat,depth) tuple of a location of interest\nγ: TMI grid\n\nOutput\n\nδ: weights on a 3D tracer field grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.linearindexXYZ-Tuple{Any}","page":"Home","title":"TMI.linearindexXYZ","text":"function linearindexXYZ(file)\nRead and assemble the grid coordinates.\n\nArguments\n\nwet: 3D mask for wet points\n\nOutput\n\nR: array of linear indices, but not a LinearIndices type\n\n\n\n\n\n","category":"method"},{"location":"#TMI.nearestneighbor","page":"Home","title":"TMI.nearestneighbor","text":"function nearestneighbor(loc,γ)\nreturn the Cartesian index and linear index \nof the nearest N neighbors\n\nArguments\n\nloc: 3-tuple of lon,lat,depth location\nγ: TMI.grid\n\nOutput\n\nInn: Cartesian indices of nearest neighbor\n\n#- Rnn: linear indices of nearest neighbor, Removed from code\n\n\n\n\n\n","category":"function"},{"location":"#TMI.nearestneighbormask","page":"Home","title":"TMI.nearestneighbormask","text":"function nearestneighbormask\nMake a 3D tracer field that is 1 at location \nof nearest neighbor, 0 elsewhere\n\nArguments\n\nloc: location in a 3-tuple (lon,lat,depth)\nγ: TMI.grid\n\nOutput\n\nδ: nearest neighbor mask 3D field\n\n\n\n\n\n","category":"function"},{"location":"#TMI.plotextent-Tuple{Any, Any}","page":"Home","title":"TMI.plotextent","text":"function plotextent\nGenerate image showing user-specified ROI\n\nArguments\n\nlatbox: in format [latstart, latstop]\nlonbox: in format [lonstart, lonstop]\n\n\n\n\n\n","category":"method"},{"location":"#TMI.readtracer-Tuple{Any, Any}","page":"Home","title":"TMI.readtracer","text":"function readtracer(file,tracername)\nRead and assemble the water-mass matrix.\n\nArguments\n\nfile: TMI NetCDF file name\ntracername: name of tracer\n\nOutput\n\nc: 3D tracer field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.regeneratedphosphate-Tuple{Any}","page":"Home","title":"TMI.regeneratedphosphate","text":"function regeneratedphosphate(TMIversion)\nRegenerated (i.e., accumulated, remineralized) phosphate\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nPO₄ᴿ: regenerated phosphate\nγ: TMI grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.sample_observations-Tuple{Any, Any, Any}","page":"Home","title":"TMI.sample_observations","text":"function sample_observations(TMIversion,variable,locs)\nSynthetic observations that are a contaminated version of real observations\nThis version: observations with random (uniform) spatial sampling\n\nArguments\n\nTMIversion::String: version of TMI water-mass/circulation model\nvariable::String: variable name to use as template\nN: number of observations\n\nOutput\n\ny: contaminated observations on 3D grid\nW⁻: appropriate weighting (inverse covariance) matrix for these observations,\nytrue: uncontaminated observations, 3D field\nlocs: 3-tuples of locations for observations\nwis: weighted indices for interpolation to locs sites\n\n\n\n\n\n","category":"method"},{"location":"#TMI.sample_observations-Tuple{Any, Any}","page":"Home","title":"TMI.sample_observations","text":"function sample_observations(TMIversion,variable)\nSynthetic observations that are a contaminated version of real observations\nThis version: gridded observations\n\nArguments\n\nTMIversion::String: version of TMI water-mass/circulation model\nvariable::String: variable name to use as template\n\nOutput\n\ny: contaminated observations on 3D grid\nW⁻: appropriate weighting (inverse covariance) matrix for these observations,\nθtrue: real observations, 3D field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.section-Tuple{Any, Any, Any}","page":"Home","title":"TMI.section","text":"function section\nView latitude-depth slice of field\n\nArguments\n\nc: 3d tracer field\nlon: longitude of section\nγ: TMI.grid\n\nOutput\n\ncsection: 2d slice of field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.steady_inversion-Union{Tuple{T}, Tuple{Vector{T}, Any, Array{T, 3}, BitArray{3}}} where T<:Real","page":"Home","title":"TMI.steady_inversion","text":"function steady_inversion(u,Alu,d,γ.wet)\ninvert for a steady-state tracer distribution\n\nArguments\n\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\nd: model constraints\nwet: BitArray ocean mask\n\nOutput\n\nc: steady-state tracer distribution\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfaceindex-Tuple{Any}","page":"Home","title":"TMI.surfaceindex","text":"function surfaceindex(I) \nGet the vector-index where depth level == 1 and it is ocean.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfaceorigin-Tuple{Any, Any}","page":"Home","title":"TMI.surfaceorigin","text":"function surfaceorigin(TMIversion,loc)\n Find the surface origin of water for some interior box \n This is equivalent to solving a sensitivity problem:\n The mass fraction at a location `loc` of interest is \n`c[loc] = δᵀ c`, where `δ` samples the location of the global mass-fraction variable, c.\nThen the sensitivity of `c[loc]` is: d(c[loc])/d(d) = A⁻ᵀ δ.\nThe derivative is solved using the constraint: Ac = d.\nThe sensitivity is exactly the mass fraction originating from each source.      \nThis problem is mathematically similar to determining how the ocean is filled.\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nloc: location (lon,lat,depth) of location of interest\n\nOutput\n\norigin: surface map of fraction of source water for a given location\nγ: TMI grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfacepatch-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, TMI.grid}} where T<:Real","page":"Home","title":"TMI.surfacepatch","text":"function surfacepatch\nMake a surface boundary condition\nwith a rectangular patch\n\nArguments\n\nlonbox: longitudes of box edges\nlatbox: latitudes of box edges\nγ: TMI.grid\n\nOutput\n\nd: vector that describes surface patch\n\n\n\n\n\n","category":"method"},{"location":"#TMI.tracerinit","page":"Home","title":"TMI.tracerinit","text":"function tracerinit(wet,ltype=Float64)\n  initialize tracer field on TMI grid\nperhaps better to have a tracer struct and constructor\n\nArguments\n\nwet::BitArray mask of ocean points\n\nOutput\n\nd:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"function"},{"location":"#TMI.trackpathways-Tuple{Any, Any, Any}","page":"Home","title":"TMI.trackpathways","text":"function trackpathways(TMIversion,latbox,lonbox)\nTrack the pathways of a user-defined water mass.\n Steps: (a) define the water mass by a rectangular surface patch dyed with passive tracer concentration of         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.\n See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nlatbox: min and max latitude of box\nlonbox: min and max longitude of box\n\nOutput\n\nc: fraction of water from surface source\nγ: TMI grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.updatelinearindex-Tuple{Any, Any, Any}","page":"Home","title":"TMI.updatelinearindex","text":"function updatelinearindex(izyx,Izyx,R)\nLinear index translated from z,y,x to x,y,z accounting\n\nArguments\n\nizyx: index of interest in z,y,x accounting\nIzyx: wet Cartesian Index for z,y,x\nR: Linear indices for x,y,z \n\nOutput\n\nixyz: index of interest in x,y,z accounting\n\n\n\n\n\n","category":"method"},{"location":"#TMI.vec2fld-Tuple{Vector{Float64}, Vector{CartesianIndex{3}}}","page":"Home","title":"TMI.vec2fld","text":"function vec2fld\nTransfer a vector to a 3D field with accounting for ocean bathymetry\n\nArguments\n\nvector: field in vector form (no land points)\nI: cartesian indices of ocean points\n\nOutput\n\nfield: field in 3d form including land points (NaN)\n\n\n\n\n\n","category":"method"},{"location":"#TMI.volumefilled-Tuple{Any}","page":"Home","title":"TMI.volumefilled","text":"function volumefilled(TMIversion)\nFind the ocean volume that has originated from each surface box.\n This is equivalent to solving a sensitivity problem:\n The total volume is V = vᵀ c , where v is the volume of each box \n and c is the fraction of volume from a given source which\n satisfies the equation A c = d.                     \n Next, dV/d(d) = A⁻ᵀ v, and dV/d(d) is exactly the volume originating from each source.\n\n See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011.\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nvolume: global ocean volume filled by a surface region\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassmatrix-Tuple{Any}","page":"Home","title":"TMI.watermassmatrix","text":"function watermassmatrix(file)\nRead and assemble the water-mass matrix.\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\nA: water-mass matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassmatrixXYZ-Tuple{Any, Any}","page":"Home","title":"TMI.watermassmatrixXYZ","text":"    function watermassmatrixXYZ(file,R)\nRead and assemble the water-mass matrix from MATLAB.\nTransfer to updated x,y,z version\n\nArguments\n\nfile: TMI NetCDF file name\nγ: TMI grid\n\nOutput\n\nA: water-mass matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassmatrixZYX-Tuple{Any}","page":"Home","title":"TMI.watermassmatrixZYX","text":"function watermassmatrixZYX(file)\nRead and assemble the water-mass matrix.\nLegacy version from MATLAB.\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\nA: water-mass matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.wetlocation-Tuple{Any}","page":"Home","title":"TMI.wetlocation","text":"function wetlocation(γ)\nGet (lon,lat,depth) tuples of wet locations.\nAllow a location to be wet if at least one out of 8 nearby gridpoints is wet.\nCertainly \"wet\" gridpoints could be defined more strictly.\n\nArguments\n\nγ: TMI.grid\n\nOutput\n\nloc: lon,lat,depth tuple\n\n\n\n\n\n","category":"method"}]
}