using Interpolations

# Error analysis
compute_rmse(field1, field2, mask, γ) = begin
    diff_tracer = replace((field1 - field2).tracer, NaN => 0.0)
    diff = diff_tracer[:, :, mask][γ.wet[:, :, mask]]
    sqrt(mean(diff.^2))
end

compute_mae(field1, field2, mask, γ) = begin
    diff_tracer = replace((field1 - field2).tracer, NaN => 0.0)
    diff = diff_tracer[:, :, mask][γ.wet[:, :, mask]]
    mean(abs.(diff))
end

compute_mean_temp(field) = mean(field)
compute_surface_temp(field, γ) = begin
    bc = getsurfaceboundary(field).tracer
    ca = cellarea(γ).tracer
    valid = .!isnan.(bc)
    sum(ca[valid] .* bc[valid]) / sum(ca[valid])
end

nan_field(field, γ) = Field(NaN * ones(size(γ.wet)), γ, field.name, field.longname, field.units)
copy_field_attributes(field, ref, γ) = Field(field.tracer, γ, ref.name, ref.longname, ref.units)
rename_field(field, γ, name) = Field(field.tracer, γ, name, field.longname, field.units)
get_z_faces(Δz) = cumsum([0, Δz...])
get_z_centers(z_faces) = (z_faces[2:end] .+ z_faces[1:end-1]) ./ 2

# Utility functions
rescale_vertical(z, (zmin_new, zmax_new)) = begin
    zmin, zmax = extrema(z)
    zmin_new == zmax_new && return fill(zmin_new, size(z))
    a = (zmax_new - zmin_new) / (zmax - zmin)
    return a .* z .+ (zmin_new - a * zmin)
end

function interpolate_with_nans(x, y, x_new; method=:linear, extrapolate_val=NaN)
    valid = .!isnan.(y)
    sum(valid) == 0 && return fill(NaN, length(x_new))
    sum(valid) == 1 && return fill(y[findfirst(valid)], length(x_new))
    length(x_new) == 1 && return fill(y[findfirst(valid)], length(x_new))

    x_v, y_v = x[valid], y[valid]
    perm = sortperm(x_v)
    x_v, y_v = x_v[perm], y_v[perm]

    extrap = extrapolate_val == :flat ? Flat() : extrapolate_val
    itp = LinearInterpolation(x_v, y_v, extrapolation_bc=extrap)
    return itp.(x_new)
end
