import Pkg; Pkg.activate("..")

# using LinearAlgebra
using TMI
# using SparseArrays
using NCDatasets

# get list of glacial inversions
glacial_versions = filter(v->occursin("LGM",v), versionlist())

for TMIversion in glacial_versions

    Amodern, Alu, γmodern, TMIfile, L, B = config(TMIversion, compute_lu = false);

    TMIfile_orig = TMIfile[1:end-3]*"_nosealeveldrop.nc"
    TMIfile_new = TMIfile[1:end-3]*"_sealeveldrop.nc"
    
    ds = NCDataset(TMIfile_orig)

    kglacial = findfirst(==(125.0),γmodern.depth)
    b_surface = ones(3, kglacial, γmodern, :bc_surface, "Surface", "nondim")

    # update grid to be consistent with boundary conditions
    #γglacial = Grid(b_surface, γmodern)
    γglacial = deepcopy(Grid(b_surface, γmodern)) # deepcopy necessary but don't know why

    # update water-mass matrix to be consistent with boundary points and grid 
    Anewboundary = watermassmatrix(Amodern, γglacial)

    # make some wet points into dry
    nx, ny = size(γglacial.wet[:,:,1])
    for k in 1:kglacial - 1
        γglacial.wet[:,:,k] = falses(nx,ny)
    end
 
    # trim water-mass matrix
    Aglacial = TMI.matrix_modern2glacial(Anewboundary, γmodern, γglacial)
    
    # update, save tracer fields
    # fnames = TMI.fieldnames()

    # update tracer fields
    for (kk,vv) in TMI.standardize_fieldnames()
        #if any(occursin.(kk,keys(ds)))    # kk in varnames || kk in xvarnames #haskey(vars,kk)
        if kk in keys(ds)    # kk in varnames || kk in xvarnames #haskey(vars,kk)
            fglacial = readfield(TMIfile_orig, kk, γglacial) # automatially updates to new grid
            writefield(TMIfile_new, fglacial)
        end
    end
    
    # for n in fnames
    #     if any(occursin.(n,keys(ds)))    

    #     end
    # end

    # up date and rewrite source
    qglacial = readsource(TMIfile_orig, "qPO₄", γglacial) # automatially updates to new grid
    TMI.writesource(TMIfile_new, qglacial)

    TMI.grid2nc(TMIversion*"_sealeveldrop", γglacial)
    
    # save matrix
    TMI.watermassmatrix2nc(TMIversion*"_sealeveldrop", Aglacial)

    # save a subset of optimization parameters
    # I'm off the hook: NetCDF has not saved these properly -- go back and see mat file
    # varname = "J" 
    # J = ncread(TMIfile_orig,"J")
    # iteratts = Dict("longname" => "iteration number")
    # Jatts =  Dict("longname" => "cost function value", "units" => "[]")
    # nccreate(TMIfile_new,varname,"iter",1:length(J),iteratts,atts=Jatts)
    
    # overwrite original file with sea-level drop version
    cp(TMIfile_new,TMIfile, force= true)
    
end


