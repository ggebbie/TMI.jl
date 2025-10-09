import Pkg; Pkg.activate(".")

using TMI
using Revise
using LinearAlgebra, Statistics, SparseArrays
using Test
using Interpolations

function rescale_thicknesses(Δz_old, H_new)
    # Calculate original total thickness
    H_old = sum(Δz_old)
    
    # Calculate scaling factor
    α = H_new / H_old
    
    # Scale each thickness
    Δz_new = α .* Δz_old
    
    return Δz_new
end

function rescale_vertical(z_old, new_bounds)
    z_min_old = minimum(z_old)
    z_max_old = maximum(z_old)
    z_min_new, z_max_new = new_bounds
    if z_min_new == z_max_new
        @warn "New bounds are equal ($z_min_new, $z_max_new). All values will be set to $z_min_new."
        return fill(z_min_new, size(z_old))
    end
    a = (z_max_new - z_min_new) / (z_max_old - z_min_old)
    b = z_min_new - a * z_min_old
    z_new = a .* z_old .+ b
    return z_new
end

# Include the interpolate function
function interpolate_with_nans(x, y, x_new; method=:linear, extrapolate_val=NaN)
    valid_idx = .!isnan.(y)
    n_valid = sum(valid_idx)
    
    if n_valid == 0
        @warn "No valid (non-NaN) points to interpolate from. Returning NaNs."
        return fill(NaN, length(x_new))
    end
    
    if n_valid == 1
        @warn "Only 1 valid point. Using constant interpolation."
        return fill(y[findfirst(valid_idx)], length(x_new))
    end
    
    x_valid = x[valid_idx]
    y_valid = y[valid_idx]
    
    sort_idx = sortperm(x_valid)
    x_valid = x_valid[sort_idx]
    y_valid = y_valid[sort_idx]
    
    # Handle the :flat symbol to use Flat() extrapolation
    if extrapolate_val == :flat
        extrapolate_val = Flat()
    end
    
    if method == :linear
        itp = LinearInterpolation(x_valid, y_valid, extrapolation_bc=extrapolate_val)
    elseif method == :cubic
        if n_valid < 4
            @warn "Cubic interpolation requires at least 4 points. Using linear instead."
            itp = LinearInterpolation(x_valid, y_valid, extrapolation_bc=extrapolate_val)
        else
            itp = CubicSplineInterpolation(x_valid, y_valid, extrapolation_bc=extrapolate_val)
        end
    elseif method == :nearest
        itp = interpolate((x_valid,), y_valid, Gridded(Constant()))
        itp = extrapolate(itp, extrapolate_val)
    else
        error("Unknown method: $method. Use :linear, :cubic, or :nearest")
    end
    
    y_new = itp.(x_new)
    
    return y_new
end

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
θ_modern_on_compressed = deepcopy(θ_compressed)

imax, jmax, kmax = size(θ_modern.tracer)
for idx = 1:imax, jdx = 1:jmax
    reference_column = θ_modern.tracer[idx, jdx, :]
    ref_wet = .!isnan.(reference_column)
    if any(.!isnan.(reference_column))

        compressed_column = θ_compressed.tracer[idx, jdx, :]
        wet_points = .!isnan.(compressed_column)
        z_centers_target = z_centers_modern[wet_points]

        if sum(wet_points) > 0
            target_bounds = extrema(z_centers_target)

            z_centers_reference_rescaled = rescale_vertical(z_centers_modern, target_bounds)
            # println(z_centers_reference_rescaled)
            reference_column_on_target = interpolate_with_nans(z_centers_reference_rescaled[ref_wet], reference_column[ref_wet], 
                                                        z_centers_target, method=:linear, extrapolate_val=Flat())

            new_compressed_column = 1 .* compressed_column
            new_compressed_column[wet_points] .= reference_column_on_target

            if sum(isnan.(new_compressed_column[wet_points] .- compressed_column[wet_points])) == 0.0
                if (new_compressed_column[wet_points][1] - reference_column[ref_wet][1]) == 0.0
                    θ_modern_on_compressed.tracer[idx, jdx, :] .= new_compressed_column
                else 
                    error("more NaNs introduced!")
                end
            else
                error("more NaNs introduced!")
            end

        end
    end
end

firstwet = findfirst(sum(γ_dropped.wet, dims = [1, 2])[:] .> 0)

using Plots

σepth = zero(γ_dropped.depth) #setup uncertainties in K  
σepth[1000 .< γ_dropped.depth] .= 0.05
σepth[γ_dropped.depth .<= 1000] .= 5.
σepth[firstwet] = 0.5
σepth[1:(firstwet-1)] .= NaN

LGM_theta_σ = get_nan_3dfield(θ_modern_on_compressed, γ_dropped) #setup uncertaity Field 
for ii in γ_dropped.I
    LGM_theta_σ.tracer[ii] = σepth[ii[3]]
end

W⁻ = (1/sum(γ_dropped.wet)) .* Diagonal(1 ./LGM_theta_σ.tracer[γ_dropped.wet].^2)
# W⁻ = Diagonal(1 ./LGM_theta_σ.tracer[γ_dropped.wet].^2)

b = zerosurfaceboundary(γ_dropped) #initial adjustment  is zero 
u = getsurfaceboundary(θ_modern_on_compressed) #first guess field 
Alu_comp = lu(A_dropped)

out, f, fg, fg! = steadyclimatology(Alu_comp,b,u,θ_modern_on_compressed,
                                    W⁻,γ_dropped; iterations = 500, method = :blackbox) #constrained line search

ũ = out.minimizer

out.res

θ̃_true = θ_modern_on_compressed 
θ̃_approx =  steadyinversion(Alu_comp, unvec(u, ũ), γ_dropped)

θ̃_approx = copy_attributes_3dfield(θ̃_approx, θ̃_true, γ_dropped)
θ̃_true = rename_3dfield(θ̃_true, γ_dropped, :θ_reference)

mean(θ̃_approx)
mean(θ̃_true)
mean(θ_modern)

TMIversion_compressed =  TMIversion*"_compressed.nc"
TMIfile_compressed = TMI.pkgdatadir("TMI_" * TMIversion_compressed)

plot(ũ - vec(u))


contourf(getsurfaceboundary(θ̃_true).tracer - getsurfaceboundary(θ_modern).tracer, clims = (-3, 3))
contourf(getsurfaceboundary(θ̃_approx).tracer - getsurfaceboundary(θ_modern).tracer, clims = (-3, 3))

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
mean(abs.(diff.^2))


depth_mask = γ_dropped.depth .== γ_dropped.depth[firstwet]
diff = θ̃_diff[:, :, depth_mask][γ_dropped.wet[:, :, depth_mask]]
sqrt(mean(diff.^2))
mean(abs.(diff.^2))

# mean(θ̃_approx.tracer[:, :, firstwet][γ_dropped.wet[:, :, firstwet]])
# mean(θ̃_true.tracer[:, :, firstwet][γ_dropped.wet[:, :, firstwet]])


# using Plots
# plot(γ_dropped.depth, θ̃_approx.tracer[50, 22, :])
# plot!(γ_dropped.depth, θ̃_true.tracer[50, 22, :])
