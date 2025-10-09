import Pkg; Pkg.activate(".")

using TMI
using Revise
using LinearAlgebra, Statistics, SparseArrays
using Test
using Interpolations
using Plots

#### HELPER FUNCTIONS #####
get_nan_3dfield(field, γ) = Field(NaN * ones(size(γ.wet)), γ, field.name, field.longname, field.units)
copy_attributes_3dfield(field, reference_field, γ) = Field(field.tracer, γ, reference_field.name, reference_field.longname, reference_field.units)
rename_3dfield(field, γ, name::Symbol) = Field(field.tracer, γ, name, field.longname, field.units)

get_z_faces(layer_thickness) = cumsum([0, layer_thickness...])
get_z_centers(z_faces) = (z_faces[2:end] .+ z_faces[1:end-1]) ./ 2.

#### HELPER FUNCTIONS #####

### BEGIN SCRIPT ####
TMIversion = "LGM_90x45x33_G14" #load in G14 that is on a modern ocean grid

drop_height = 120.0 #select the height that you would like to compress the ocean by
compression_lower_boundary = Inf #select the lower boundary of compression

A_modern, Alu, γ_modern, TMIfile_modern, L, B = config(TMIversion * "_nosealeveldrop");
A_dropped, Alu, γ_dropped, TMIfile_dropped, L, B = config(TMIversion);

Δz_modern = layerthickness(γ_modern)
z_faces_modern = get_z_faces(Δz_modern)
z_centers_modern = get_z_centers(z_faces_modern)
num_z_levels_modern = length(z_centers_modern)

θ_compressed = readfield(TMIfile_dropped, "θ", γ_dropped) #read temperature field on compressed grid
θ_modern = readfield(TMIfile_modern, "θ", γ_modern) #read temperature field on modern grid

firstwet = findfirst(sum(γ_dropped.wet, dims = [1, 2])[:] .> 0)
θ_modern_on_compressed = deepcopy(θ_compressed)
θ_modern_on_compressed.tracer[:, :, firstwet] .= (0 .* getsurfaceboundary(θ_modern_on_compressed).tracer) .+ getsurfaceboundary(θ_modern).tracer


σepth = zero(γ_dropped.depth) #setup uncertainties in K  
σepth[1000 .< γ_dropped.depth] .= 0.01
σepth[γ_dropped.depth .<= 1000] .= 20.
σepth[firstwet] = 0.01
σepth[1:(firstwet-1)] .= NaN

LGM_theta_σ = get_nan_3dfield(θ_modern_on_compressed, γ_dropped) #setup uncertaity Field 
for ii in γ_dropped.I
    LGM_theta_σ.tracer[ii] = σepth[ii[3]]
end

W⁻ = (1/sum(γ_dropped.wet)) .* Diagonal(1 ./LGM_theta_σ.tracer[γ_dropped.wet].^2)

b = zerosurfaceboundary(γ_dropped) #initial adjustment  is zero 
u = getsurfaceboundary(θ_modern_on_compressed) #first guess field 
Alu_comp = lu(A_dropped)

out, f, fg, fg! = steadyclimatology(Alu_comp,b,u,θ_modern_on_compressed,
                                    W⁻,γ_dropped; iterations = 1000, method = :gradient) #constrained line search

ũ = out.minimizer
plot(ũ - vec(u))                         


θ̃_true = θ_modern_on_compressed 
θ̃_approx =  steadyinversion(Alu_comp, unvec(u, ũ), γ_dropped)

θ̃_approx = copy_attributes_3dfield(θ̃_approx, θ̃_true, γ_dropped)
θ̃_true = rename_3dfield(θ̃_true, γ_dropped, :θ_reference)

mean(θ̃_approx)
mean(θ_modern)

area = vec(cellarea(γ_dropped))
area_notdrop = vec(cellarea(γ_modern))

sum(vec(getsurfaceboundary(θ̃_approx)) .* area) / sum(area)
sum(vec(getsurfaceboundary(θ_modern_on_compressed)) .* area) / sum(area)
sum(vec(getsurfaceboundary(θ_compressed)) .* area) / sum(area)
sum(vec(getsurfaceboundary(θ_modern)) .* area_notdrop) / sum(area_notdrop)

mean(vec(getsurfaceboundary(θ_modern)))

plot(ũ - vec(u))

heatmap(getsurfaceboundary(θ_modern_on_compressed).tracer' - getsurfaceboundary(θ_modern).tracer', clims = (-3, 3))
heatmap(getsurfaceboundary(θ̃_approx).tracer' - getsurfaceboundary(θ_modern).tracer', clims = (-3, 3))

heatmap(getsurfaceboundary(θ̃_approx).tracer')
heatmap(getsurfaceboundary(θ̃_true).tracer')


TMIversion_compressed =  TMIversion*"_compressed.nc"
TMIfile_compressed = TMI.pkgdatadir("TMI_" * TMIversion_compressed)

isfile(TMIfile_compressed) && rm(TMIfile_compressed)

writefield(TMIfile_compressed, θ̃_approx)
writefield(TMIfile_compressed, θ̃_true)

TMI.grid2nc(TMIversion*"_compressed", γ_dropped)
TMI.watermassmatrix2nc(TMIversion*"_compressed", A_dropped)


v = cellvolume(γ_dropped).tracer

θ̃_diff = replace((θ̃_true - θ̃_approx).tracer, NaN => 0.0)

depth_mask = 1000 .<= γ_dropped.depth .<= Inf
diff = θ̃_diff[:, :, depth_mask][γ_dropped.wet[:, :, depth_mask]]
sqrt(mean(diff.^2))
mean(abs.(diff))


depth_mask = 200 .<= γ_dropped.depth .<= 1000
diff = θ̃_diff[:, :, depth_mask][γ_dropped.wet[:, :, depth_mask]]
sqrt(mean(diff.^2))
mean(abs.(diff))

depth_mask = γ_dropped.depth .== γ_dropped.depth[firstwet]
diff = θ̃_diff[:, :, depth_mask][γ_dropped.wet[:, :, depth_mask]]
sqrt(mean(diff.^2))
mean(abs.(diff))

mean(θ̃_approx.tracer[:, :, firstwet][γ_dropped.wet[:, :, firstwet]])
mean(θ̃_true.tracer[:, :, firstwet][γ_dropped.wet[:, :, firstwet]])


using Plots
plot(γ_dropped.depth, θ̃_approx.tracer[50, 22, :])
plot!(γ_dropped.depth, θ̃_true.tracer[50, 22, :])
