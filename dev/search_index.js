var documenterSearchIndex = {"docs":
[{"location":"#","page":"Home","title":"Home","text":"CurrentModule = TMI","category":"page"},{"location":"#TMI-1","page":"Home","title":"TMI","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Documentation for TMI.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [TMI]","category":"page"},{"location":"#Base.:*-Union{Tuple{T}, Tuple{Any, TMI.BoundaryCondition{T}}} where T<:Real","page":"Home","title":"Base.:*","text":"`function *(C,d::BoundaryCondition)::BoundaryCondition`\nDefine scalar or matrix multiplication for BoundaryCondition`s\n\n\n\n\n\n","category":"method"},{"location":"#Base.:*-Union{Tuple{T}, Tuple{Any, TMI.Field{T}}} where T<:Real","page":"Home","title":"Base.:*","text":"`function *(C,d::Field)::Field`\nDefine scalar or matrix multiplication for fields\n\n\n\n\n\n","category":"method"},{"location":"#Base.:*-Union{Tuple{T}, Tuple{TMI.Field{T}, TMI.Field{T}}} where T<:Real","page":"Home","title":"Base.:*","text":"`function *(c::Field,d::Field)::Field`\nField by field multiplication is element-by-element.\n\n\n\n\n\n","category":"method"},{"location":"#Base.:+-Union{Tuple{T}, Tuple{TMI.BoundaryCondition{T}, TMI.BoundaryCondition{T}}} where T<:Real","page":"Home","title":"Base.:+","text":"`function +(c::BoundaryCondition,d::BoundaryCondition)::BoundaryCondition`\nDefine addition for Fields\n\n\n\n\n\n","category":"method"},{"location":"#Base.:+-Union{Tuple{T}, Tuple{TMI.Field{T}, TMI.Field{T}}} where T<:Real","page":"Home","title":"Base.:+","text":"`function +(c::Field,d::Field)::Field`\nDefine addition for Fields\n\n\n\n\n\n","category":"method"},{"location":"#Base.:--Union{Tuple{T}, Tuple{TMI.Field{T}, TMI.Field{T}}} where T<:Real","page":"Home","title":"Base.:-","text":"`function -(c::Field,d::Field)::Field`\nDefine addition for Fields\n\n\n\n\n\n","category":"method"},{"location":"#Base.ones-Tuple{Int64, Int64, BitArray{3}}","page":"Home","title":"Base.ones","text":"Initialize boundary condition with ones\n\n\n\n\n\n","category":"method"},{"location":"#Base.zeros","page":"Home","title":"Base.zeros","text":"function zeros(wet,ltype=Float64)\ninitialize tracer field on TMI grid\nThis version will give an array\n\nArguments\n\nwet::BitArray mask of ocean points\nltype:: optional type argument, default=Float64\n\nOutput\n\nd:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"function"},{"location":"#Base.zeros-Tuple{Int64, Int64, BitArray{3}}","page":"Home","title":"Base.zeros","text":"Initialize boundary condition with zeroes\n\n\n\n\n\n","category":"method"},{"location":"#Base.zeros-Tuple{TMI.Grid}","page":"Home","title":"Base.zeros","text":"function zeros(γ::Grid)\n  initialize tracer field on TMI grid\n  using a Field struct and constructor\n\nArguments\n\nγ::TMI.Grid\n\nOutput\n\nd::Field,  3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"method"},{"location":"#LinearAlgebra.dot-Union{Tuple{T}, Tuple{TMI.Field{T}, TMI.Field{T}}} where T<:Real","page":"Home","title":"LinearAlgebra.dot","text":"`function *(c::Field,d::Field)::Field`\nField by field multiplication is element-by-element.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.boundarymatrix-Tuple{Any, Any}","page":"Home","title":"TMI.boundarymatrix","text":"    function boundarymatrix(file,γ)\nRead and assemble the boundary matrix from MATLAB.\nTransfer to updated x,y,z version\n\nArguments\n\nfile: TMI MATLAB file name\nγ: TMI grid\n\nOutput\n\nB: boundary condition matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.cartesianindex-Tuple{BitArray{3}}","page":"Home","title":"TMI.cartesianindex","text":"function cartesianindex(wet)\nRead and assemble the grid coordinates\naccording to a 3D tracer in x,y,z order\n\nArguments\n\nwet: BitArray logical mask for wet points\n\nOutput\n\nI: 3D Cartesian indices\n\n\n\n\n\n","category":"method"},{"location":"#TMI.cartesianindex-Tuple{String}","page":"Home","title":"TMI.cartesianindex","text":"function cartesianindex(file)\nRead and assemble the grid coordinates\naccording to the legacy MATLAB code (z,y,x order).\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\nI: TMI Cartesian index for wet points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.cellarea-Tuple{Any}","page":"Home","title":"TMI.cellarea","text":"Horizontal area of grid cell\n\n\n\n\n\n","category":"method"},{"location":"#TMI.cellvolume-Tuple{Any}","page":"Home","title":"TMI.cellvolume","text":"Volume of each grid cell.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.ces_ncwrite-Tuple{Any, Any, Any}","page":"Home","title":"TMI.ces_ncwrite","text":"function ces_ncwrite(γ,time,sol_array)\nWrite .nc file output for commonerasim.jl\n\nArguments\n\nγ: \ntime: vector of time values \nsol_array: solution array in form time x lat x lon x depth - must match γ + time \n\nOutput\n\nsaves .nc file titled \"ces_output.nc\" in data array \n\n\n\n\n\n","category":"method"},{"location":"#TMI.circulationmatrix-Tuple{Any, Any, Any}","page":"Home","title":"TMI.circulationmatrix","text":"function circulationmatrix(file,A,γ)\nRead and assemble the circulation matrix from the efficient storage of A and F₀ variables.\n\nArguments\n\nfile: TMI MATLAB file name\nA: TMI water-mass matrix\nγ: TMI grid\n\nOutput\n\nL: circulation matrix in xyz format\n\n\n\n\n\n","category":"method"},{"location":"#TMI.circulationmatrix-Tuple{Any, Any}","page":"Home","title":"TMI.circulationmatrix","text":"function circulationmatrix(file,γ)\nRead and assemble the circulation matrix from MATLAB.\nTransfer to updated x,y,z version\n\nArguments\n\nfile: TMI MATLAB file name\nγ: TMI grid\n\nOutput\n\nL: circulation matrix in xyz format\n\n\n\n\n\n","category":"method"},{"location":"#TMI.config2nc-NTuple{5, Any}","page":"Home","title":"TMI.config2nc","text":"Save TMI configuration to NetCDF format for non-proprietary access\n\n\n\n\n\n","category":"method"},{"location":"#TMI.config_from_mat-Tuple{Any}","page":"Home","title":"TMI.config_from_mat","text":"Configure TMI environment from original MATLAB output\n\n\n\n\n\n","category":"method"},{"location":"#TMI.config_from_nc-Tuple{Any}","page":"Home","title":"TMI.config_from_nc","text":"function config_from_nc(TMIversion)\nConfigure TMI environment from NetCDF input file.\n\nArguments\n\nTMIversion: TMI version for water-mass/circulation model\n\nOutput\n\nA: TMI steady-state water-mass matrix\nAlu: LU decomposition of A\nγ: TMI grid properties\nTMIfile: TMI file name\n\n\n\n\n\n","category":"method"},{"location":"#TMI.control2state-Union{Tuple{T}, Tuple{Matrix{T}, Any}} where T<:Real","page":"Home","title":"TMI.control2state","text":"function control2state(tracer2D,γ)\nturn 2D surface field into 3D field with zeroes below surface\n\nArguments\n\ntracer2D:: 2D surface tracer field\nwet::BitArray mask of ocean points\n\nOutput\n\ntracer3D:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction!-Union{Tuple{T}, Tuple{Any, Any, Vector{T}, Any, TMI.BoundaryCondition{T}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, Any, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.costfunction!","text":"function costfunction!(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,Q⁻,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\nIssue: couldn't figure out how to nest with costfunction_obs!\nIssue: why are wis and locs both needed? `gobserve` function\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction-Union{Tuple{T}, Tuple{Vector{T}, Any, TMI.BoundaryCondition{T}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, Any, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.costfunction","text":"function costfunction(J,gJ,uvec,Alu,b,y,Wⁱ,wis,Q⁻,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\nIssue #1: couldn't figure out how to nest with costfunction_obs!\nIssue #2: Update for BoundaryCondition types\n\nArguments\n\nuvec: controls, vector format\nAlu: LU decomposition of water-mass matrix\nb: boundary condition\ny: pointwise observations\nWⁱ: inverse of W weighting matrix for observations\nwis: weights for interpolation (data sampling, E)\nlocs: data locations (lon,lat,depth)\nQ⁻: weights for control vector\nγ: grid\n\nOutput\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs!-Union{Tuple{T}, Tuple{Any, Any, Vector{T}, Any, Array{T, 3}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs!","text":"function costfunction_obs!(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,locs,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ndfld: model constraints\ny: pointwise observations\nWⁱ: inverse of W weighting matrix for observations\nwis: weights for interpolation (data sampling, E)\nlocs: data locations (lon,lat,depth)\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs!-Union{Tuple{T}, Tuple{Any, Any, Vector{T}, Any, TMI.BoundaryCondition{T}, TMI.Field{T}, LinearAlgebra.Diagonal{T, Vector{T}}, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs!","text":"function costfunction_obs!(J,gJ,u::BoundaryCondition{T},Alu,b::BoundaryCondition{T},y::Field{T},Wⁱ::Diagonal{T, Vector{T}}) where T <: Real\n\nsquared model-data misfit for gridded data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\nu: controls, field format\nAlu: LU decomposition of water-mass matrix\nb: boundary conditions\ny: observations on grid\nWⁱ: inverse of W weighting matrix for observations\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs-Union{Tuple{T}, Tuple{Vector{T}, Any, Array{T, 3}, Vector{T}, LinearAlgebra.Diagonal{T, Vector{T}}, Any, Any, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs","text":"function costfunction_obs(u,Alu,dfld,yfld,Wⁱ,wis,locs,γ)\nsquared model-data misfit for pointwise data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nu: controls, vector format\nAlu: LU decomposition of water-mass matrix\ny: pointwise observations\nd: model constraints\nWⁱ: inverse of W weighting matrix for observations\nwis: weights for interpolation \nlocs: data locations\nγ: grid\n\nOutput\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\n\n\n\n\n\n","category":"method"},{"location":"#TMI.costfunction_obs-Union{Tuple{T}, Tuple{Vector{T}, Any, TMI.BoundaryCondition{T}, TMI.Field{T}, LinearAlgebra.Diagonal{T, Vector{T}}, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.costfunction_obs","text":"function costfunction_obs(uvec::Vector{T},Alu,b::BoundaryCondition{T},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where T <: Real\n\nsquared model-data misfit for gridded data\ncontrols are a vector input for Optim.jl\n\nArguments\n\nJ: cost function of sum of squared misfits\ngJ: derivative of cost function wrt to controls\nu: controls, field format\nAlu: LU decomposition of water-mass matrix\nb: boundary conditions\ny: observations on grid\nWⁱ: inverse of W weighting matrix for observations\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.depthindex-Tuple{Any}","page":"Home","title":"TMI.depthindex","text":"function depthindex(I) \nGet the k-index (depth level) from the Cartesian index\n\n\n\n\n\n","category":"method"},{"location":"#TMI.fld2vec-Union{Tuple{T}, Tuple{Array{T, 3}, Vector{CartesianIndex{3}}}} where T<:Real","page":"Home","title":"TMI.fld2vec","text":"function fld2vec\nTransfer 3D field with accounting for ocean bathymetry to a vector without land points.\nThis is done more easily with a BitArray mask, i.e., vector = field[mask].\nThis function may be removed in the future.\n\nArguments\n\nfield: field in 3d form including land points (NaN)\nI: cartesian indices of ocean points\n\nOutput\n\nvector: field in vector form (no land points)\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gridprops-Tuple{Any}","page":"Home","title":"TMI.gridprops","text":"function gridprops(file)\nRead and assemble the grid properties.\n\nArguments\n\nfile: TMI NetCDF file name\n\nOutput\n\ngrid: TMI grid coordinates\n\n\n\n\n\n","category":"method"},{"location":"#TMI.horizontaldistance-Tuple{Any, TMI.Grid}","page":"Home","title":"TMI.horizontaldistance","text":"function horizontaldistance(loc,γ)\nreturn the Cartesian index and linear index \nof the nearest N neighbors\n\nArguments\n\nloc: 3-tuple of lon,lat,depth location\nγ: TMI.grid\n\nOutput\n\nhordist: horizontal distance to nearest tracer grid points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.interpindex-Tuple{Any, Any}","page":"Home","title":"TMI.interpindex","text":"function interpindex(loc,γ)     Weights for linear interpolation.     The derivative of linear interpolation is needed in sensitivity studies.     ReverseDiff.jl could find this quantity automatically.     Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.\n\nArguments\n\nc: a temporary tracer field, would be nice to make it unnecessary\nloc: (lon,lat,depth) tuple of a location of interest\nγ: TMI grid\n\nOutput\n\nδ: weights on a 3D tracer field grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.interpweights-Tuple{Any, Any}","page":"Home","title":"TMI.interpweights","text":"function interpweights(loc,γ)     Weights for linear interpolation.     The derivative of linear interpolation is needed in sensitivity studies.     ReverseDiff.jl could find this quantity automatically.     Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.\n\nArguments\n\nloc: (lon,lat,depth) tuple of a location of interest\nγ: TMI grid\n\nOutput\n\nδ: weights on a 3D tracer field grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.latindex-Tuple{Any}","page":"Home","title":"TMI.latindex","text":"function latindex(I) \nGet the j-index (latitude index) from the Cartesian index\n\n\n\n\n\n","category":"method"},{"location":"#TMI.linearindex-Tuple{Any}","page":"Home","title":"TMI.linearindex","text":"function linearindex(wet)\nRead and assemble the grid coordinates.\n\nArguments\n\nwet: 3D mask for wet points\n\nOutput\n\nR: array of linear indices, but not a LinearIndices type\n\n\n\n\n\n","category":"method"},{"location":"#TMI.lonindex-Tuple{Any}","page":"Home","title":"TMI.lonindex","text":"function lonindex(I) \nGet the i-index (lon index) from the Cartesian index\n\n\n\n\n\n","category":"method"},{"location":"#TMI.matrix_zyx2xyz-Tuple{Any, Any, Any}","page":"Home","title":"TMI.matrix_zyx2xyz","text":"    function matrix_zyx2xyz(TMIfile,Azyx,γ)\n\nTransfer zyx format water-mass matrix A to xyz format\n\nArguments\n\nAzyx: water-mass matrix in zyx format\nγ: TMI grid\n\nOutput\n\nAxyz: water-mass matrix in xyz format\n\n\n\n\n\n","category":"method"},{"location":"#TMI.nearestneighbor","page":"Home","title":"TMI.nearestneighbor","text":"function nearestneighbor(loc,γ)\nreturn the Cartesian index and linear index \nof the nearest N neighbors\n\nArguments\n\nloc: 3-tuple of lon,lat,depth location\nγ: TMI.grid\n\nOutput\n\nInn: Cartesian indices of nearest neighbor\n\n#- Rnn: linear indices of nearest neighbor, Removed from code\n\n\n\n\n\n","category":"function"},{"location":"#TMI.nearestneighbormask","page":"Home","title":"TMI.nearestneighbormask","text":"function nearestneighbormask\nMake a 3D tracer field that is 1 at location \nof nearest neighbor, 0 elsewhere\n\nArguments\n\nloc: location in a 3-tuple (lon,lat,depth)\nγ: TMI.grid\n\nOutput\n\nδ: nearest neighbor mask 3D field\n\n\n\n\n\n","category":"function"},{"location":"#TMI.observe-Union{Tuple{T}, Tuple{TMI.Field{T}, Array{Tuple{Interpolations.WeightedAdjIndex{2, T}, Interpolations.WeightedAdjIndex{2, T}, Interpolations.WeightedAdjIndex{2, T}}, 1}, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.observe","text":"function observe\nTake a observation at location given by weights wis\n\n\n\n\n\n","category":"method"},{"location":"#TMI.oxygen-Tuple{Any, Any}","page":"Home","title":"TMI.oxygen","text":"Reconstruct dissolved oxygen (that doesn't exist in TMI product) by assuming same oxygen saturation fraction as modern\n\n\n\n\n\n","category":"method"},{"location":"#TMI.planviewplot-Union{Tuple{T}, Tuple{TMI.BoundaryCondition{T}, Any, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.planviewplot","text":"function planviewplot\nPlot of plan view (lon-lat) in ocean\n\nArguments\n\nfield::BoundaryCondition, 3d filed of values to be plotted\ndepth: depth of plan view\nlims: contour levels\nγ::Grid, needed for lat, lon but not in BoundaryCondition! (could refactor)\ntitlelabel: optional title label\n\n\n\n\n\n","category":"method"},{"location":"#TMI.planviewplot-Union{Tuple{T}, Tuple{TMI.Field{T}, Any, Any}} where T<:Real","page":"Home","title":"TMI.planviewplot","text":"function planviewplot\nPlot of plan view (lon-lat) in ocean\n\nArguments\n\nfield::Field, 3d filed of values to be plotted\ndepth: depth of plan view\nlims: contour levels\ntitlelabel: optional title label\n\n\n\n\n\n","category":"method"},{"location":"#TMI.plotextent-Tuple{Any, Any}","page":"Home","title":"TMI.plotextent","text":"function plotextent\nGenerate image showing user-specified ROI\n\nArguments\n\nlatbox: in format [latstart, latstop]\nlonbox: in format [lonstart, lonstop]\n\n\n\n\n\n","category":"method"},{"location":"#TMI.readfield-Tuple{Any, Any, TMI.Grid}","page":"Home","title":"TMI.readfield","text":"function readfield(file,tracername,γ)\nRead a tracer field from NetCDF but return it \nas a Field.\n\nArguments\n\nfile: TMI NetCDF file name\ntracername: name of tracer\nγ::Grid, TMI grid specification\n\nOutput\n\nc::Field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.readtracer-Tuple{Any, Any}","page":"Home","title":"TMI.readtracer","text":"function readtracer(file,tracername)\nRead a tracer field from NetCDF.\n\nArguments\n\nfile: TMI NetCDF file name\ntracername: name of tracer\n\nOutput\n\nc: 3D tracer field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.regeneratedphosphate-Tuple{Any, Any, Any}","page":"Home","title":"TMI.regeneratedphosphate","text":"function regeneratedphosphate(TMIversion,Alu,γ)\nRegenerated (i.e., accumulated, remineralized) phosphate\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nAlu: LU decomposition of water-mass matrix A\nγ: TMI grid\n\nOutput\n\nPO₄ᴿ: regenerated phosphate\n\n\n\n\n\n","category":"method"},{"location":"#TMI.section-Union{Tuple{T}, Tuple{TMI.Field{T}, Any}} where T<:Real","page":"Home","title":"TMI.section","text":"function section\nView latitude-depth slice of field\n\nArguments\n\nc::Field, 3D tracer field plus meta data\nlon: longitude of section\n\nOutput\n\ncsection: 2d slice of field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.sectionplot-Union{Tuple{T}, Tuple{TMI.Field{T}, Any, Any}} where T<:Real","page":"Home","title":"TMI.sectionplot","text":"function sectionplot\nPlot of section (lat-depth) in ocean\n\nArguments\n\nfield::Field, 3d filed of values to be plotted\nlon: longitude of section\nlims: contour levels\ntitlelabel: optional title labeln\n\n\n\n\n\n","category":"method"},{"location":"#TMI.setboundarycondition!-Union{Tuple{T}, Tuple{TMI.Field{T}, TMI.BoundaryCondition{T}}} where T<:Real","page":"Home","title":"TMI.setboundarycondition!","text":"function setboundarycondition!(d::Field,b::BoundaryCondition)\napply boundary condition to the equation constraints\n\nArguments\n\nd::Field, equation constraints (i.e., right hand side)\nb::BoundaryCondition\n\n\n\n\n\n","category":"method"},{"location":"#TMI.setsource!-Union{Tuple{T}, Tuple{TMI.Field{T}, TMI.Field{T}}, Tuple{TMI.Field{T}, TMI.Field{T}, Any}} where T<:Real","page":"Home","title":"TMI.setsource!","text":"function setsource!(d::Field,q::Field,r::Number)\napply interior source q to the equation constraints d\n\nArguments\n\nd::Field, equation constraints (i.e., right hand side)\nq::Field, interior source\nr::Number, default = 1.0, stoichiometric ratio\n\n\n\n\n\n","category":"method"},{"location":"#TMI.sparsedatamap-Union{Tuple{T}, Tuple{Vector{T}, Any, TMI.BoundaryCondition{T}, Vector{T}, Any, Array{Tuple{Interpolations.WeightedAdjIndex{2, T}, Interpolations.WeightedAdjIndex{2, T}, Interpolations.WeightedAdjIndex{2, T}}, 1}, Any, Any, TMI.Grid}, Tuple{Vector{T}, Any, TMI.BoundaryCondition{T}, Vector{T}, Any, Array{Tuple{Interpolations.WeightedAdjIndex{2, T}, Interpolations.WeightedAdjIndex{2, T}, Interpolations.WeightedAdjIndex{2, T}}, 1}, Any, Any, TMI.Grid, Any}} where T<:Real","page":"Home","title":"TMI.sparsedatamap","text":"function sparsedatamap(u₀::Vector{T},Alu,b::BoundaryCondition{T},y::Vector{T},W⁻,wis,locs,Q⁻,γ::Grid;iterations=10) where T <: Real\n\n Find the distribution of a tracer given:\n (a) the pathways described by A or its LU decomposition Alu,\n (b) first-guess boundary conditions and interior sources given by d₀,\n (c) perturbations to the surface boundary condition u₀\nthat best fits observations, y,\naccording to the cost function,\nJ = (ỹ - y)ᵀ W⁻¹ (ỹ - y)\nsubject to Aỹ = d₀ + Γ u₀.                 \nW⁻ is a (sparse) weighting matrix.\nSee Supplementary Section 2, Gebbie & Huybers 2011.\n\nArguments\n\nu₀:\nAlu:\nb: first guess of boundary conditions and interior sources\ny: observations on 3D grid\nW⁻: weighting matrix best chosen as inverse error covariance matrix\nfg!: compute cost function and gradient in place\nγ: grid\n\n\n\n\n\n","category":"method"},{"location":"#TMI.steadyclimatology-Tuple{Any, Any, Any}","page":"Home","title":"TMI.steadyclimatology","text":"function steadyclimatology(u₀,fg!,iterations)      Find the distribution of a tracer given:      (a) the pathways described by A or its LU decomposition Alu,      (b) first-guess boundary conditions and interior sources given by d₀,      (c) perturbations to the surface boundary condition u₀     that best fits observations, y,     according to the cost function,     J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)     subject to Aỹ = d₀ + Γ u₀.                      W⁻ is a (sparse) weighting matrix.     See Supplementary Section 2, Gebbie & Huybers 2011.\n\nArguments\n\nu₀:\nfg!: compute cost function and gradient in place\niterations: number of optimization iterations\n\n\n\n\n\n","category":"method"},{"location":"#TMI.steadyinversion-Union{Tuple{T}, Tuple{Any, TMI.BoundaryCondition{T}, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.steadyinversion","text":"function steadyinversion(Alu,b;q=nothing,r=1.0)\ninvert for a steady-state tracer distribution\n\nArguments\n\nAlu: LU decomposition of water-mass matrix\nb: boundary condition\nγ::Grid\n\nOptional Arguments\n\nq: interior sources/sinks of phosphate\nr: stochiometric ratio of tracer:phosphate\n\nOutput\n\nc::Field, steady-state tracer distribution\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surface_oxygensaturation-Tuple{Any}","page":"Home","title":"TMI.surface_oxygensaturation","text":"Surface oxygen saturation value and fraction of saturation value in field \n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfacecontrol2field!-Union{Tuple{T}, Tuple{Array{T, 3}, Vector{T}, Any}} where T<:Real","page":"Home","title":"TMI.surfacecontrol2field!","text":"function surfacecontrol2field!(c,u,γ)\nAdd surface control vector to existing 3D field\n\nArguments\n\nc:: state field, 3d tracer field with NaN on dry points, modified by function\nusfc:: surface control vector\nwet::BitArray mask of ocean points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfacecontrol2field!-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, Any}} where T<:Real","page":"Home","title":"TMI.surfacecontrol2field!","text":"function surfacecontrol2field!(c,u,γ)\nAdd surface control vector to tracer vector\n\nArguments\n\nc:: state field, 3d tracer field with NaN on dry points, modified by function\nu:: surface control vector\nwet::BitArray mask of ocean points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfacecontrol2field-Union{Tuple{T}, Tuple{Vector{T}, Any}} where T<:Real","page":"Home","title":"TMI.surfacecontrol2field","text":"function surfacecontrol2field(usfc,γ.wet)\nturn surface control vector into 3D field with zeroes below surface\n\nArguments\n\nusfc:: surface control vector\nwet::BitArray mask of ocean points\n\nOutput\n\ntracer3D:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfaceindex-Tuple{Any}","page":"Home","title":"TMI.surfaceindex","text":"function surfaceindex(I) \nGet the vector-index where depth level == 1 and it is ocean.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfaceorigin-Tuple{Any, Any, TMI.Grid}","page":"Home","title":"TMI.surfaceorigin","text":"function surfaceorigin(TMIversion,loc)\n Find the surface origin of water for some interior box \n This is equivalent to solving a sensitivity problem:\n The mass fraction at a location `loc` of interest is \n`c[loc] = δᵀ c`, where `δ` samples the location of the global mass-fraction variable, c.\nThen the sensitivity of `c[loc]` is: d(c[loc])/d(d) = A⁻ᵀ δ.\nThe derivative is solved using the constraint: Ac = d.\nThe sensitivity is exactly the mass fraction originating from each source.      \nThis problem is mathematically similar to determining how the ocean is filled.\n\nArguments\n\nloc: location (lon,lat,depth) of location of interest\nAlu: LU decomposition of water-mass matrix A\nγ: TMI grid\n\nOutput\n\norigin: surface map of fraction of source water for a given location, log10 of effective depth, in terms of a BoundaryCondition\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfacepatch-Tuple{Any, Any, TMI.Grid}","page":"Home","title":"TMI.surfacepatch","text":"function surfacepatch\nMake a surface boundary condition\nwith a rectangular patch\n\nArguments\n\nlonbox: longitudes of box edges\nlatbox: latitudes of box edges\nγ: TMI.grid\n\nOutput\n\nd: vector that describes surface patch\n\n\n\n\n\n","category":"method"},{"location":"#TMI.synthetic_observations-NTuple{4, Any}","page":"Home","title":"TMI.synthetic_observations","text":"function synthetic_observations(TMIversion,variable,locs)\nSynthetic observations that are a contaminated version of real observations\nThis version: observations with random (uniform) spatial sampling\n\nArguments\n\nTMIversion::String: version of TMI water-mass/circulation model\nvariable::String: variable name to use as template\nN: number of observations\n\nOutput\n\ny: contaminated observations on 3D grid\nW⁻: appropriate weighting (inverse covariance) matrix for these observations,\nytrue: uncontaminated observations, 3D field\nlocs: 3-tuples of locations for observations\nwis: weighted indices for interpolation to locs sites\n\n\n\n\n\n","category":"method"},{"location":"#TMI.synthetic_observations-Tuple{Any, Any, Any}","page":"Home","title":"TMI.synthetic_observations","text":"function synthetic_observations(TMIversion,variable)\nSynthetic observations that are a contaminated version of real observations\nThis version: gridded observations\n\nArguments\n\nTMIversion::String: version of TMI water-mass/circulation model\nvariable::String: variable name to use as template\n\nOutput\n\ny: contaminated observations on 3D grid\nW⁻: appropriate weighting (inverse covariance) matrix for these observations,\nθtrue: real observations, 3D field\n\n\n\n\n\n","category":"method"},{"location":"#TMI.tracerinit-Tuple{Any, Any, Any}","page":"Home","title":"TMI.tracerinit","text":"function tracerinit(wet,vec,I)\n      initialize tracer field on TMI grid\n    perhaps better to have a tracer struct and constructor\n\nArguments\n\nwet:: BitArray mask of ocean points\nvec:: vector of values at wet points\nI:: Cartesian Index for vector\n\nOutput\n\nfield:: 3d tracer field with NaN on dry points\n\n\n\n\n\n","category":"method"},{"location":"#TMI.trackpathways-NTuple{4, Any}","page":"Home","title":"TMI.trackpathways","text":"function trackpathways(TMIversion,latbox,lonbox)\nTrack the pathways of a user-defined water mass.\n Steps: (a) define the water mass by a rectangular surface patch dyed with passive tracer concentration of         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.\n See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nlatbox: min and max latitude of box\nlonbox: min and max longitude of box\nγ: TMI grid\n\nOutput\n\nc: fraction of water from surface source\n\n\n\n\n\n","category":"method"},{"location":"#TMI.updatelinearindex-Tuple{Any, Any, Any}","page":"Home","title":"TMI.updatelinearindex","text":"function updatelinearindex(izyx,Izyx,R)\nLinear index translated from z,y,x to x,y,z accounting\n\nArguments\n\nizyx: index of interest in z,y,x accounting\nIzyx: wet Cartesian Index for z,y,x\nR: Linear indices for x,y,z \n\nOutput\n\nixyz: index of interest in x,y,z accounting\n\n\n\n\n\n","category":"method"},{"location":"#TMI.varying!-NTuple{4, Any}","page":"Home","title":"TMI.varying!","text":"function varying!(du, u, p, t)\nODE function for varying boundary cond\nSets up dc/dt = L*C + B*f to be solved\n\nArguments\n\ndu: dc/dt (must have this name for DifferentialEquations.jl to work\nu: C, what we are solving for \np: parameters for diffeq - must hold specified vars  \nt: time we are solving for (automatically determined by DE.jl)\n\nOutput\n\ndu: numerical value of LC+Bf, vector of size 74064 for 4°\n\n\n\n\n\n","category":"method"},{"location":"#TMI.vec2fld-Union{Tuple{T}, Tuple{Vector{T}, Vector{CartesianIndex{3}}}} where T<:Real","page":"Home","title":"TMI.vec2fld","text":"function vec2fld\nTransfer a vector to a 3D field with accounting for ocean bathymetry\n\nArguments\n\nvector: field in vector form (no land points)\nI: cartesian indices of ocean points\n\nOutput\n\nfield: field in 3d form including land points (NaN)\n\n\n\n\n\n","category":"method"},{"location":"#TMI.volumefilled-Tuple{Any, Any, Any}","page":"Home","title":"TMI.volumefilled","text":"function volumefilled(TMIversion)\nFind the ocean volume that has originated from each surface box.\n This is equivalent to solving a sensitivity problem:\n The total volume is V = vᵀ c , where v is the volume of each box \n and c is the fraction of volume from a given source which\n satisfies the equation A c = d.                     \n Next, dV/d(d) = A⁻ᵀ v, and dV/d(d) is exactly the volume originating from each source.\n\n See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011.\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nAlu: LU decomposition of water-mass matrix A\nγ: TMI.grid\n\nOutput\n\nvolume: log10 of global ocean volume filled by a surface region, exists at surface, therefore given BoundaryCondition type\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassdistribution-NTuple{4, Any}","page":"Home","title":"TMI.watermassdistribution","text":"function watermassdistribution(TMIversion,latbox,lonbox)\nTrack the pathways of a user-defined water mass.\n Steps: (a) define the water mass by an oceanographically-relevant surface patch dyed with passive tracer concentration of one\n     (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.\n See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\nAlu: LU decomposition of water-mass matrix A\nregion: name of pre-defined surface region\nγ: TMI grid\n\nOutput\n\ng: water-mass fraction\n\n\n\n\n\n","category":"method"},{"location":"#TMI.watermassmatrix-Tuple{Any}","page":"Home","title":"TMI.watermassmatrix","text":"function watermassmatrix(file)\nRead and assemble the water-mass matrix.\n\nArguments\n\nfile: TMI NetCDF or MATLAB file name\n\nOutput\n\nA: water-mass matrix\n\n\n\n\n\n","category":"method"},{"location":"#TMI.wetlocation-Tuple{Any}","page":"Home","title":"TMI.wetlocation","text":"function wetlocation(γ)\nGet (lon,lat,depth) tuples of wet locations.\nAllow a location to be wet if at least one out of 8 nearby gridpoints is wet.\nCertainly \"wet\" gridpoints could be defined more strictly.\n\nArguments\n\nγ: TMI.grid\n\nOutput\n\nloc: lon,lat,depth \n\n\n\n\n\n","category":"method"},{"location":"#TMI.BoundaryCondition","page":"Home","title":"TMI.BoundaryCondition","text":"struct BoundaryCondition\n\na plane defined at `dim=dimval`\nCan array have other element types?\nAre indices needed?\n\n\n\n\n\n","category":"type"},{"location":"#TMI.Field","page":"Home","title":"TMI.Field","text":"struct Tracer\n\nThis structure permits the grid to be \nautomatically passed to functions with\nthe tracer field.\n\nThis structure assumes the Tracer type to be \nthree-dimensional.\n\n\n\n\n\n","category":"type"},{"location":"#Base.:\\-Union{Tuple{T}, Tuple{Any, TMI.Field{T}}} where T<:Real","page":"Home","title":"Base.:\\","text":"`function \\(A,d::Field)::Field`\nDefine left division for Fields\nNeed two slashes to prevent invalid escape\n\n\n\n\n\n","category":"method"},{"location":"#TMI.E","page":"Home","title":"TMI.E","text":"function E \nE anonymously calls field2obs\n\n\n\n\n\n","category":"function"},{"location":"#TMI.boundarymatrix2nc-Tuple{Any, Any}","page":"Home","title":"TMI.boundarymatrix2nc","text":"Save boundary matrix for transient model to NetCDF file\n\n\n\n\n\n","category":"method"},{"location":"#TMI.circulationmatrix2nc-Tuple{Any, Any, Any}","page":"Home","title":"TMI.circulationmatrix2nc","text":"Save circulation matrix L to NetCDF file.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.eastindex-Tuple{Any}","page":"Home","title":"TMI.eastindex","text":"function eastindex(I) \nGet the vector index on the northern open boundary\n\n\n\n\n\n","category":"method"},{"location":"#TMI.fieldsatts-Tuple{}","page":"Home","title":"TMI.fieldsatts","text":"All variable names and attributes. Useful for writing NetCDF files.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.getboundarycondition-NTuple{4, Any}","page":"Home","title":"TMI.getboundarycondition","text":"Get boundary condition by extracting from 3D tracer\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gobserve-Union{Tuple{T}, Tuple{Vector{T}, TMI.Field{T}, Any}} where T<:Real","page":"Home","title":"TMI.gobserve","text":"function gobserve(gy::Vector{T},c::Field{T},wis,γ) where T <: Real\n\nADJOINT Take a observation at location given by weights wis\nArguments not symmetric with `observe` due to splat operator\n\n\n\n\n\n","category":"method"},{"location":"#TMI.grid2nc-Tuple{Any, Any}","page":"Home","title":"TMI.grid2nc","text":"Put grid properties (Cartesian index) into NetCDF file\n\n\n\n\n\n","category":"method"},{"location":"#TMI.griddicts-Tuple{Any}","page":"Home","title":"TMI.griddicts","text":"Save grid dictionaries of attributes for writing to NetCDF file\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gsetboundarycondition-Union{Tuple{T}, Tuple{TMI.Field{T}, TMI.BoundaryCondition{T}}} where T<:Real","page":"Home","title":"TMI.gsetboundarycondition","text":"function gsetboundarycondition(gd::Field{T},b::BoundaryCondition{T}) where T<: Real\n\nADJOINT: apply boundary condition to the equation constraints\n\nArguments\n\nd::Field, equation constraints (i.e., right hand side)\nb::BoundaryCondition\n\n\n\n\n\n","category":"method"},{"location":"#TMI.gsteadyinversion-Union{Tuple{T}, Tuple{Any, Any, TMI.BoundaryCondition{T}, TMI.Grid}} where T<:Real","page":"Home","title":"TMI.gsteadyinversion","text":"function steadyinversion(Alu,b;q=nothing,r=1.0)\ninvert for a steady-state tracer distribution\n\nArguments\n\nAlu: LU decomposition of water-mass matrix\nb: boundary condition\nγ::Grid\n\nOptional Arguments\n\nq: interior sources/sinks of phosphate\nr: stochiometric ratio of tracer:phosphate\n\nOutput\n\nc::Field, steady-state tracer distribution\n\n\n\n\n\n","category":"method"},{"location":"#TMI.matfields2nc-Tuple{Any, Any}","page":"Home","title":"TMI.matfields2nc","text":"Read 3D fields from mat file and save to NetCDF file.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.maturl-Tuple{Any}","page":"Home","title":"TMI.maturl","text":"function maturl(TMIversion)\nFind *mat file here.\nplaceholder function to give location (URL) of Google Drive input\nin the future, consider a struct or Dict that describes all TMI versions.\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nurl: location (URL) for download\n\n\n\n\n\n","category":"method"},{"location":"#TMI.ncurl-Tuple{Any}","page":"Home","title":"TMI.ncurl","text":"function ncurl(TMIversion)\nplaceholder function to give location (URL) of NetCDF Google Drive input\nin the future, consider a struct or Dict that describes all TMI versions.\n\nArguments\n\nTMIversion: version of TMI water-mass/circulation model\n\nOutput\n\nurl: location (URL) for download\n\n\n\n\n\n","category":"method"},{"location":"#TMI.northindex-Tuple{Any}","page":"Home","title":"TMI.northindex","text":"function northindex(I) \nGet the vector index on the northern open boundary\n\n\n\n\n\n","category":"method"},{"location":"#TMI.optim2nc-Tuple{Any}","page":"Home","title":"TMI.optim2nc","text":"Save optimization parameters to NetCDF file)\n\nFuture considerations: split into 2 functions\n\nread from mat\nsave to nc\n\n\n\n\n\n","category":"method"},{"location":"#TMI.regions2nc-Tuple{Any, Any}","page":"Home","title":"TMI.regions2nc","text":"Read vectors from mat file, translate to 3D,  and save surface field to NetCDF file.\n\n\n\n\n\n","category":"method"},{"location":"#TMI.southindex-Tuple{Any}","page":"Home","title":"TMI.southindex","text":"function southindex(I) \nGet the vector-index on the southern open boundary\n\n\n\n\n\n","category":"method"},{"location":"#TMI.surfaceregion-Tuple{String, String, TMI.Grid}","page":"Home","title":"TMI.surfaceregion","text":" function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition\n\nRead an oceanographically-relevant surface region from NetCDF file. (Also could be read from mat file.)\nReturn a BoundaryCondition\n\n\n\n\n\n","category":"method"},{"location":"#TMI.westindex-Tuple{Any}","page":"Home","title":"TMI.westindex","text":"function westindex(I) \nGet the vector index on the western open boundary\n\n\n\n\n\n","category":"method"},{"location":"#TMI.Γsfc","page":"Home","title":"TMI.Γsfc","text":"function Γsfc \nΓsfc anonymously calls surfacecontrol2field\n\n\n\n\n\n","category":"function"},{"location":"#TMI.Γsfc!","page":"Home","title":"TMI.Γsfc!","text":"function Γsfc! \nΓsfc! anonymously calls surfacecontrol2field!\n\n\n\n\n\n","category":"function"}]
}
