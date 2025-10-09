

# =============================================================================
# EXAMPLE: Using both functions together
# =============================================================================

println("Example: Rescaling and Interpolating Vertical Grid Data\n")

# Original data: depths in meters with some NaN values
z_old_m = [-1.0, 0.0, 25.0, 50.0, 75.0, 100.0, 110.0]
temperature = [NaN, 20.0, NaN, 15.0, 12.0, 10.0, NaN]  # Temperature data with NaN

# Step 2: Rescale depths from meters to feet (0-100m → 0-328ft)
z_old_ft = rescale_vertical(z_old_m, (0, 328))

# Step 3: Create new uniform depth grid in feet
z_new_ft = 0:32.8:328  # Every ~33 feet

# Step 4: Interpolate temperature values to new grid

temperature_new = interpolate_with_nans(z_old_ft, temperature, z_old_ft, 
                                        method=:linear, extrapolate_val=NaN)