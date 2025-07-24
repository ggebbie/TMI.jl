import Pkg; Pkg.activate(".")

using TMI
using Revise
using LinearAlgebra, Statistics, SparseArrays
using Test

get_nan_3dfield(x, γ) = Field(NaN * ones(size(γ.wet)),γ,x.name,x.longname,x.units)
get_z_faces(Δz) = cumsum([0, Δz...])
get_z_centers(zf) = (zf[2:end] .+ zf[1:end-1]) ./ 2.
function get_errors(θ̃_approx, θ̃_true, γ_compressed, depth_mask)
    rms(x, y) = sqrt(mean((x .- y).^2))
    abserr(x, y) = mean(abs.(x .- y))
    wet_points = γ_compressed.wet .* reshape(depth_mask, 1, 1, 33)
    model = θ̃_approx.tracer[wet_points]
    data = θ̃_true.tracer[wet_points]

    return (rms_err = rms(model, data), abs_err = abserr(model, data))
end

TMIversion = "LGM_90x45x33_G14" #load in G14 that is on a modern ocean grid
Amodern, Alu, γ_modern, TMIfile, L, B = config(TMIversion * "_nosealeveldrop");

drop_height = 120.0 #select the height that you would like to compress the ocean by
compression_lower_boundary = 1000.0 #select the lower boundary of compression


Δz_modern = layerthickness(γ_modern)
zf_modern = get_z_faces(Δz_modern)
zc_modern = get_z_centers(zf_modern)
nz_modern = length(zc_modern)

##calculate new z-centers based on compression parameters##
compression_boundary_index = Base.findfirst(zc_modern .>= compression_lower_boundary) - 1

H_upper = sum(Δz_modern[1:compression_boundary_index])
compression_scaling = (H_upper - drop_height) / H_upper

Δz_compressed = copy(Δz_modern)
Δz_compressed[1:compression_boundary_index] .*= compression_scaling

zf_compressed = get_z_faces(Δz_compressed)
zc_compressed = get_z_centers(zf_compressed)

##generate a new wet mask by using modern ocean depth as a criteria##
deptho_modern = sum((γ_modern.wet .* reshape(Δz_modern, 1, 1, nz_modern)), dims = 3)[:, :, 1] #calculate the ocean depth on modern grid
drop_height_mask = (deptho_modern .>= drop_height) #find layers that cannot by compressed 
wet_compressed = (drop_height_mask .* γ_modern.wet) #generate new wet mask 

##generate a compressed grid##
b_surface_compressed = deepcopy(ones(3, 1, γ_modern, :bc_surface, "Surface", "nondim"))
b_surface_compressed.tracer .= replace(x -> x == 1 ? 1.0 : NaN, 1.0 .* wet_compressed[:, :, 1]) #update the boundary conditions 
b_surface_compressed.wet .= copy(wet_compressed[:, :, 1]) #update the wet mask for boundary condition 

γ_compressed = deepcopy(Grid(b_surface_compressed, γ_modern)) # deepcopy necessary but don't know why
γ_compressed.depth .= copy(zc_compressed) #update z centers
γ_compressed.wet .= copy(wet_compressed) #update wet mask

## Modify A matrix to conform to new geometry 
rows_modern = γ_modern.R
rows_compressed_in_modern = rows_modern[γ_compressed.wet]
rows_modern = rows_modern[γ_modern.wet]

A_comp = copy(Amodern)[rows_compressed_in_modern, rows_compressed_in_modern] #subset A matrix with ocean points in compressed grid
rows = γ_compressed.R[γ_compressed.wet]
new_land_rows = rows[sum(A_comp[rows, :], dims = 2) .< -1e-14] #determine where mass fractions do not sum to 1

for landrow in new_land_rows #update mass fractions in A matrix if previously contained a shallow point
    tmprow = copy(A_comp[landrow, :])
    tmprow[landrow] = 0.0 #zero out the diagonal entry
    tmprow ./= sum(tmprow) #renormalize mass fractions
    tmprow[landrow] = -1.0
    A_comp[landrow, :] = tmprow
end

@test minimum(γ_compressed.R[γ_compressed.wet]) == 1
@test maximum(γ_compressed.R) == sum(γ_compressed.wet)
@test isapprox(maximum(sum(A_comp,dims=2)),1.0)
@test minimum(sum(A_comp,dims=2))> -1e-14

Alu_comp = lu(A_comp)

## Update tracer fields
θtrue_on_compressed = readfield(TMIfile,"θ",γ_compressed) #read the old temperature field on the new grid

σepth = zero(γ_compressed.depth) #setup uncertainties in K  
σepth[1000 .< γ_compressed.depth] .= 0.01
σepth[γ_compressed.depth .<= 1000] .= 0.5
σepth[1] = 5

LGM_theta_σ = get_nan_3dfield(θtrue_on_compressed, γ_compressed)
for ii in γ_compressed.I
    LGM_theta_σ.tracer[ii] = σepth[ii[3]]
end
W⁻ = (1/sum(γ_compressed.wet)) .* Diagonal(1 ./LGM_theta_σ.tracer[γ_compressed.wet].^2)

b = zerosurfaceboundary(γ_compressed) #initial adjustment  is zero 
u = getsurfaceboundary(θtrue_on_compressed) #first guess field 

out, f, fg, fg! = steadyclimatology(Alu_comp,b,u,θtrue_on_compressed,W⁻,γ_compressed; iterations = 100)

ũ = out.minimizer
θ̃_approx = steadyinversion(Alu_comp, unvec(u, ũ), γ_compressed)
θ̃_true = θtrue_on_compressed 

TMIversion_compressed =  TMIversion*"_compressed.nc"
TMIfile_compressed = TMI.pkgdatadir("TMI_" * TMIversion_compressed)

isfile(TMIfile_compressed) && rm(TMIfile_compressed)

writefield(TMIfile_compressed, θ̃_approx)
TMI.grid2nc(TMIversion*"_compressed", γ_compressed)
TMI.watermassmatrix2nc(TMIversion*"_compressed", A_comp)

compressed_depths = γ_compressed.depth
err_dict = Dict()
err_dict["surface"] = get_errors(θ̃_approx, θ̃_true, γ_compressed, compressed_depths .< 3)
err_dict["0to1000"] = get_errors(θ̃_approx, θ̃_true, γ_compressed, 3 .<  compressed_depths .<= 1000)
err_dict["1000tobottom"] = get_errors(θ̃_approx, θ̃_true, γ_compressed, 1000 .<  compressed_depths)
err_dict

minimum(ũ)
mean(θ̃_approx)
mean(θ̃_true)

# function get_v(Alu, γ)
#     area = cellarea(γ)
#     v = volumefilled("",Alu,γ)
#     v = (10 .^ v.tracer) .* area.tracer
# end

# v = get_v(A_compressedlu, γ_compressed)
# sum(v[.!isnan.(v)]) == sum(v[γ_compressed.wet[:, :, 1]])
# v = v[γ_compressed.wet[:, :, 1]]
# (sum(v) - sum(cellvolume(γ_compressed))) / sum(cellvolume(γ_compressed))

# v = get_v(Amodern, γ_modern)
# sum(v[.!isnan.(v)]) == sum(v[γ_modern.wet[:, :, 1]])
# v = v[γ_modern.wet[:, :, 1]]
# (sum(v) - sum(cellvolume(γ_modern))) / sum(cellvolume(γ_modern))


# LGM_θ = readfield(TMIfile,"θ",γ_modern)
# LGM_SST = getsurfaceboundary(LGM_θ)

# LGM_SST_dropped = deepcopy(b_surface_compressed)
# LGM_SST_dropped.tracer .= LGM_SST.tracer .* LGM_SST_dropped.tracer

# mean(LGM_SST.tracer[.!isnan.(LGM_SST.tracer)])
# mean(LGM_SST_dropped.tracer[.!isnan.(LGM_SST_dropped.tracer)])

# LGM_θ_dropped = steadyinversion(A_compressedlu, LGM_SST_dropped, γ_compressed)

# mean(readfield(TMIfile,"θ",γ_compressed))
# mean(LGM_θ_dropped)
# mean(LGM_θ)

# θtrue_on_compressed = readfield(TMIfile,"θ",γ_compressed)

# # u = LGM_SST_dropped #zerosurfaceboundary(γ_compressed)
# # b =  getsurfaceboundary(θtrue_on_compressed)



# mean(ũ)
# mean(vec(LGM_SST_dropped))

# unvec(u, ũ) - LGM_SST_dropped

# mean(steadyinversion(A_compressedlu, unvec(u, ũ), γ_compressed))
# mean(LGM_θ)
# mean(LGM_θ_dropped)

# J,∂J∂u = fg(ũ)
# J₀,∂J∂u0 = fg(vec(u))
# @test J < J₀
# @test out.minimum < J₀
# @test isapprox(J,out.minimum)

# ∇f, ∇f_finite = TMI.gradient_check(vec(u),f,fg,fg!)
# @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1



