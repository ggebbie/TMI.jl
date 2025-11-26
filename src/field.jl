"""
    struct Field

    This structure permits the grid to be 
    automatically passed to functions with
    the tracer field.

    This structure assumes the Tracer type to be 
    three-dimensional.

    tracer::AbstractArray{T,N}
    γ::Grid{A,N}
    name::Symbol
    longname::String
    units::String
"""
struct Field{T <: Real,R <: Real,N,F <: AbstractArray{T,N}}
    tracer::F
    γ::Grid{R,N}
    name::Symbol
    longname::String
    units::String
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, x::Field)
    summary(io, x); println(io)
    print(io, "Field size ")
    println(io, size(x.tracer))
    println(io, "Surface view")
    naxes = length(x.γ.axes)
    if  naxes == 3
        show(io,mime,heatmap(transpose(x.tracer[:,:,surfaceindex(x.γ)]),zlabel=x.units,title=x.longname))
    elseif naxes == 1
        show(io,mime,heatmap(transpose(x.tracer[:]),zlabel=x.units,title=x.longname))
    else 
        println("heatmap not supported for Field with axes length $naxes")
    end
end

""" 
    function zeros(γ::Grid,name=:none,longname="unknown",units="unknown")::Field

      initialize tracer field on TMI grid
      using a Field struct and constructor
# Arguments
- `γ`::TMI.Grid
# Output
- `d`::Field,  3d tracer field with NaN on dry points
"""
function zeros(γ::Grid{T},name=:none,longname="unknown",units="unknown")::Field where T <: Real

    # preallocate
    tracer = Array{T}(undef,size(γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer[γ.wet] .= zero(T)
    tracer[.!γ.wet] .= zero(T) / zero(T) # NaNs with right type
    return Field(tracer,γ,name,longname,units)
end

function planview(c::Field{T},depth)::Matrix{T} where T <: Real
 
    isec = findall(==(depth),c.γ.depth)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    cplan = dropdims(view(c.tracer,:,:,isec),dims=3)
    return cplan
end

"""
    function section
    View latitude-depth slice of field
# Arguments
- `c::Field`, 3D tracer field plus meta data
- `lon`: longitude of section
# Output
- `csection`: 2d slice of field
"""
function section(c::Field{T},lon)::Array{T,2} where T <: Real

    # handle longitudinal ambiguities
    if lon < minimum(c.γ.lon)
        lon += 360
    end
    
    isec = findall(==(lon),c.γ.lon)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    csection= dropdims(view(c.tracer,isec,:,:),dims=1)
    return csection
end

"""
    function readfield(file,tracername,γ)
    Read a tracer field from NetCDF but return it 
    as a Field.

    Use NCDatasets so that Unicode is correct

# Arguments
- `file`: TMI NetCDF file name
- `tracername`: name of tracer
- `γ::Grid`, TMI grid specification
# Output
- `c`::Field

---------------------------------------------------
    MATLAB version
    function readfield(matfile,mattracername,γ::Grid,Izyx) # for MATLAB

    read MATLAB field and transfer zyx format to xyz
"""
function readfield(file,tracername,γ::Grid{A,N}) where {A,N} 

    # The mode "r" stands for read-only. The mode "r" is the default mode and the parameter can be omitted.
    tracer, units, longname = _read3d(file,tracername)
    T = eltype(tracer)
    checkgrid!(tracer,γ.wet)
    c = Field(tracer,γ,tracerdict()[tracername],longname,units)
    return c
end
function readfield(matfile,mattracername,γ::Grid,Izyx) 
    # read MATLAB field and transfer zyx format to xyz
    matobj = matopen(matfile)
    varnames, xvarnames = matvarnames(matfile)

    if mattracername in varnames
        tvar = read(matobj,mattracername)
    elseif mattracername in xvarnames
        tvar = read(matobj,"x")[mattracername]
    end

    # put zyx vector into xyz 3D array
    tracer = tracerinit(tvar, Izyx, γ.wet)
    checkgrid!(tracer,γ.wet)

    nctracername = mat2ncfield()[mattracername]
    units = fieldsatts()[nctracername]["units"]
    longname = fieldsatts()[nctracername]["longname"]

    close(matobj)
    return Field(tracer,γ,Symbol(nctracername),longname,units)
end
readmatfield(file,mattracername,γ::Grid,Izyx = cartesianindex(file)) = readfield(file,mattracername,γ,Izyx)

writefield(file,c::Field) = write(file,c) 

# complete list of TMI field names
fieldnames() = ("Θ", "θ", "σθ",
                "S⋆", "Sₚ", "σSₚ",
                "δ¹⁸Ow", "σδ¹⁸Ow",
                "PO₄", "σPO₄",
                "NO₃", "σNO₃",
                "O₂", "σO₂",
                "δ¹³C", "σδ¹³C")

"""
    function mat2ncfield

    Rename MATLAB variables to NetCDF variables
"""
standardize_fieldnames() = Dict("CT"=>"Θ", "Θ"=>"Θ",
    "θ"=>"θ","Tobs"=>"θ","Tmod"=>"θ","Tlgm"=>"θ",
    "Terr"=>"σθ", "σθ"=>"σθ",
    "Sstar"=>"S⋆", "S⋆"=>"S⋆",
    "Sobs"=>"Sₚ","Smod"=>"Sₚ","Slgm"=>"Sₚ","Sₚ"=>"Sₚ", "Sp" => "Sₚ", 
    "Serr"=>"σSₚ", "σSₚ"=>"σSₚ",
    "O18obs"=>"δ¹⁸Ow","O18mod"=>"δ¹⁸Ow","O18lgm"=>"δ¹⁸Ow","δ¹⁸Ow"=>"δ¹⁸Ow",
    "O18err"=>"σδ¹⁸Ow", "σδ¹⁸Ow"=>"σδ¹⁸Ow",
    "Pobs"=>"PO₄","Pmod"=>"PO₄","Plgm"=>"PO₄","P"=>"PO₄", "PO₄"=>"PO₄",
    "Perr" => "σPO₄", "σPO₄"=> "σPO₄",
    "Nobs"=>"NO₃","Nmod"=>"NO₃","Nlgm"=>"NO₃","N"=>"NO₃","NO₃"=>"NO₃",
    "Nerr" => "σNO₃","σNO₃" => "σNO₃",
    "Oobs"=>"O₂","Omod"=>"O₂","Olgm"=>"O₂","O"=>"O₂","O₂"=>"O₂",
    "Oerr"=>"σO₂", "σO₂"=>"σO₂",
    "C13obs"=>"δ¹³C","C13mod"=>"δ¹³C","C13lgm"=>"δ¹³C","δ¹³C"=>"δ¹³C",
    "C13err" =>  "σδ¹³C", "σδ¹³C" =>  "σδ¹³C")

standardize_sourcenames() = Dict("dP"=>"qPO₄","q"=>"qPO₄","qPO₄"=>"qPO₄")
