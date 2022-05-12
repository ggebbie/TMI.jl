var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TMI","category":"page"},{"location":"#TMI","page":"Home","title":"TMI","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TMI.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TMI]","category":"page"},{"location":"#TMI.cartesianindexXYZ-Tuple{Any}","page":"Home","title":"TMI.cartesianindexXYZ","text":"function cartesianindexXYZ(wet)\nRead and assemble the grid coordinates\naccording to a 3D tracer in x,y,z order\n\nArguments\n\nwet: BitArray logical mask for wet points\n\nOutput\n\nI: 3D Cartesian indices\n\n\n\n\n\n","category":"method"},{"location":"#TMI.cartesianindexZYX-Tuple{Any}","page":"Home","title":"TMI.cartesianindexZYX","text":"function cartesianindexZYX(file)\nRead and assemble the grid coordinates\naccording to the legacy MATLAB code (z,y,x order).\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\ngrid: TMI grid coordinates\n\n\n\n\n\n","category":"method"},{"location":"#TMI.config-Tuple{Any}","page":"Home","title":"TMI.config","text":"function config(url,inputdir)\n\nArguments\n\nTMIversion: TMI version for water-mass/circulation model\n\nOutput\n\nA: TMI steady-state water-mass matrix\nAlu: LU decomposition of A\nγ: TMI grid properties\nTMIfile: TMI file name\n\n\n\n\n\n","category":"method"},{"location":"#TMI.depthindex-Tuple{Any}","page":"Home","title":"TMI.depthindex","text":"function depthindex(I) \n\nGet the k-index (depth level) from the Cartesian index\n\n\n\n\n\n","category":"method"},{"location":"#TMI.download-Tuple{Any, Any}","page":"Home","title":"TMI.download","text":"function download(url,inputdir)\nRead and assemble all TMI inputs.\n\nArguments\n\nurl: Google Drive location of TMI input\ninputdir: input directory location to store file\n\nOutput\n\nnone\n\n\n\n\n\n","category":"method"},{"location":"#TMI.dyeplot-NTuple{4, Any}","page":"Home","title":"TMI.dyeplot","text":"function dyeplot\nPlot of dye in ocean\n\nArguments\n\nlat: latitude arrays\ndepth: depth array\nvals: lat x depth value array\nlims: contour levels\n\n\n\n\n\n","category":"method"},{"location":"#TMI.fld2vec-Tuple{Array{Float64, 3}, Vector{CartesianIndex{3}}}","page":"Home","title":"TMI.fld2vec","text":"function fld2vec\nTransfer 3D field with accounting for ocean bathymetry to a vector without land points\n\nArguments\n\nfield: field in 3d form including land points (NaN)\nI: cartesian indices of ocean points\n\nOutput\n\nvector: field in vector form (no land points)\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gdriveurl-Tuple{Any}","page":"Home","title":"TMI.gdriveurl","text":"function gdriveurl(TMIversion)\nplaceholder function to give location (URL) of Google Drive input\nin the future, consider a struct or Dict that describes all TMI versions.\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nurl: location (URL) for download\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gridprops-Tuple{Any}","page":"Home","title":"TMI.gridprops","text":"function gridprops(file)\nRead and assemble the grid properties.\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\ngrid: TMI grid coordinates\n\n\n\n\n\n","category":"method"},{"location":"#TMI.horizontaldistance-Tuple{Any, TMI.grid}","page":"Home","title":"TMI.horizontaldistance","text":"function horizontaldistance(loc,γ)\nreturn the Cartesian index and linear index \nof the nearest N neighbors\n\nArguments\n\nloc: 3-tuple of lon,lat,depth location\nγ: TMI.grid\n\nOutput\n\nhordist: horizontal distance to nearest tracer grid points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.linearindexXYZ-Tuple{Any}","page":"Home","title":"TMI.linearindexXYZ","text":"function linearindexXYZ(file)\nRead and assemble the grid coordinates.\n\nArguments\n\nwet: 3D mask for wet points\n\nOutput\n\nR: array of linear indices, but not a LinearIndices type\n\n\n\n\n\n","category":"method"},{"location":"#TMI.misfit_gridded_data-Union{Tuple{T}, Tuple{Vector{T}, Any, Array{T, 3}, Array{T, 3}, LinearAlgebra.Diagonal{T, Vector{T}}, BitArray{3}}} where T<:Real","page":"Home","title":"TMI.misfit_gridded_data","text":"function misfit_gridded_data(u,Alu,y,d,Wⁱ,wet)\nsquared model-data misfit\ncontrols are a vector input for Optim.jl\n\nArguments\n\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ny: observations on grid\nd: model constraints\nWⁱ: inverse of W weighting matrix for observations\nwet: BitArray ocean mask\n\nOutput\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\n\n\n\n\n\n","category":"method"},{"location":"#TMI.nearestneighbor-Tuple{Any, Any}","page":"Home","title":"TMI.nearestneighbor","text":"function nearestneighbor(loc,γ)\nreturn the Cartesian index and linear index \nof the nearest N neighbors\n\nArguments\n\nloc: 3-tuple of lon,lat,depth location\nγ: TMI.grid\n\nOutput\n\nInn: Cartesian indices of nearest neighbor\nRnn: linear indices of nearest neighbor\n\n\n\n\n\n","category":"method"},{"location":"#TMI.nearestneighbormask-Tuple{Any, TMI.grid}","page":"Home","title":"TMI.nearestneighbormask","text":"function nearestneighbormask\nMake a 3D tracer field that is 1 at location \nof nearest neighbor, 0 elsewhere\n\nArguments\n\nloc: location in a 3-tuple (lon,lat,depth)\nγ: TMI.grid\n\nOutput\n\nδ: nearest neighbor mask 3D field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.plotextent-Tuple{Any, Any}","page":"Home","title":"TMI.plotextent","text":"function plotextent\nGenerate image showing user-specified ROI\n\nArguments\n\nlatbox: in format [latstart, latstop]\nlonbox: in format [lonstart, lonstop]\n\n\n\n\n\n","category":"method"},{"location":"#TMI.readtracer-Tuple{Any, Any}","page":"Home","title":"TMI.readtracer","text":"function readtracer(file,tracername)\nRead and assemble the water-mass matrix.\n\nArguments\n\nfile: TMI NetCDF file name\ntracername: name of tracer\n\nOutput\n\nc: 3D tracer field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.regeneratedphosphate-Tuple{Any}","page":"Home","title":"TMI.regeneratedphosphate","text":"function regeneratedphosphate(TMIversion)\nRegenerated (i.e., accumulated, remineralized) phosphate\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nPO₄ᴿ: regenerated phosphate\nγ: TMI grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.section-Tuple{Any, Any, Any}","page":"Home","title":"TMI.section","text":"function section\nView latitude-depth slice of field\n\nArguments\n\nc: 3d tracer field\nlon: longitude of section\nγ: TMI.grid\n\nOutput\n\ncsection: 2d slice of field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfaceindex-Tuple{Any}","page":"Home","title":"TMI.surfaceindex","text":"function surfaceindex(I) \n\nGet the vector-index where depth level == 1 and it is ocean.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfaceorigin-Tuple{Any, Any}","page":"Home","title":"TMI.surfaceorigin","text":"function surfaceorigin(TMIversion,loc)\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nloc: location (lon,lat,depth) of location of interest\n\nOutput\n\norigin: surface map of fraction of source water for a given location\nγ: TMI grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfacepatch-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, TMI.grid}} where T<:Real","page":"Home","title":"TMI.surfacepatch","text":"function surfacepatch\nMake a surface boundary condition\nwith a rectangular patch\n\nArguments\n\nlonbox: longitudes of box edges\nlatbox: latitudes of box edges\nγ: TMI.grid\n\nOutput\n\nd: vector that describes surface patch\n\n\n\n\n\n","category":"method"},{"location":"#TMI.tracerinit","page":"Home","title":"TMI.tracerinit","text":"function tracerinit(wet,ltype=Float64)\n  initialize tracer field on TMI grid\nperhaps better to have a tracer struct and constructor\n\nArguments\n\nwet::BitArray mask of ocean points\n\nOutput\n\nd:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"function"},{"location":"#TMI.trackpathways-Tuple{Any, Any, Any}","page":"Home","title":"TMI.trackpathways","text":"function trackpathways(TMIversion,latbox,lonbox)\nFraction of water originating from surface box\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nlatbox: min and max latitude of box\nlonbox: min and max longitude of box\n\nOutput\n\nc: fraction of water from surface source\nγ: TMI grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.updatelinearindex-Tuple{Any, Any, Any}","page":"Home","title":"TMI.updatelinearindex","text":"function updatelinearindex(izyx,Izyx,R)\nLinear index translated from z,y,x to x,y,z accounting\n\nArguments\n\nizyx: index of interest in z,y,x accounting\nIzyx: wet Cartesian Index for z,y,x\nR: Linear indices for x,y,z \n\nOutput\n\nixyz: index of interest in x,y,z accounting\n\n\n\n\n\n","category":"method"},{"location":"#TMI.vec2fld-Tuple{Vector{Float64}, Vector{CartesianIndex{3}}}","page":"Home","title":"TMI.vec2fld","text":"function vec2fld\nTransfer a vector to a 3D field with accounting for ocean bathymetry\n\nArguments\n\nvector: field in vector form (no land points)\nI: cartesian indices of ocean points\n\nOutput\n\nfield: field in 3d form including land points (NaN)\n\n\n\n\n\n","category":"method"},{"location":"#TMI.volumefilled-Tuple{Any}","page":"Home","title":"TMI.volumefilled","text":"function volumefilled(TMIversion)\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nvolume: global ocean volume filled by a surface region\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassmatrix-Tuple{Any}","page":"Home","title":"TMI.watermassmatrix","text":"function watermassmatrix(file)\nRead and assemble the water-mass matrix.\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\nA: water-mass matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassmatrixXYZ-Tuple{Any, Any}","page":"Home","title":"TMI.watermassmatrixXYZ","text":"    function watermassmatrixXYZ(file,R)\nRead and assemble the water-mass matrix from MATLAB.\nTransfer to updated x,y,z version\n\nArguments\n\nfile: TMI NetCDF file name\nγ: TMI grid\n\nOutput\n\nA: water-mass matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassmatrixZYX-Tuple{Any}","page":"Home","title":"TMI.watermassmatrixZYX","text":"function watermassmatrixZYX(file)\nRead and assemble the water-mass matrix.\nLegacy version from MATLAB.\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\nA: water-mass matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.Γ-Union{Tuple{T}, Tuple{Matrix{T}, Any}} where T<:Real","page":"Home","title":"TMI.Γ","text":"function Γ(tracer2D,γ)\nturn 2D surface field into 3D field with zeroes below surface\n\nArguments\n\ntracer2D:: 2D surface tracer field\nwet::BitArray mask of ocean points\n\nOutput\n\ntracer3D:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"method"}]
}