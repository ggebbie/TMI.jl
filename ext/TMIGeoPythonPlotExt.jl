module TMIGeoPythonPlotExt

using TMI
using GeoPythonPlot
import GeoPythonPlot: planviewplot

"""
function `planviewplot`

Plot of plan view (lon-lat) in ocean using cartopy and PythonPlot (matplotlib)

# Arguments
- `field::Field`, 3d filed of values to be plotted
- `depth`: depth of plan view
- `lims`: contour levels
# Optional Arguments
- `titlelabel`: optional title label
- `fname`: output file name
- `cenlon`: central longitude
"""
function GeoPythonPlot.planviewplot(c::Field{<:Real}, depth, lims;titlelabel="planview plot",fname="fname.png",cenlon=-160.0) 

    cplan = planview(c,depth)
    lon = c.γ.lon
    lat = c.γ.lat
    
    # replace Field with Array
    GeoPythonPlot.planviewplot(cplan, lon, lat, depth, lims; titlelabel=titlelabel,fname=fname,cenlon=cenlon,units=c.units) 

end

# #TMI.sectionplot()
# """
#     function sectionplot
#     Plot of section (lat-depth) in ocean
# # Arguments
# - `field::Field`, 3d filed of values to be plotted
# - `lon`: longitude of section
# - `lims`: contour levels
# - `titlelabel`: optional title labeln
# """
# function sectionplot(field::Field{<:Real}, lon, lims;titlelabel="section plot",fname="figure.png") 

#     Psection = section(field,lon)
#     cmap_seismic = get_cmap("seismic")
#     z = field.γ.depth/1000.0
    
#     #calc fignum - based on current number of figures
#     figure()
#     #cmap_seismic.set_bad(color="black") # doesn't work
#     contourf(field.γ.lat, z, Psection', lims, cmap=cmap_seismic)
#     colorbar(orientation="horizontal")
#     CS = gca().contour(field.γ.lat, z, Psection', lims,colors="k")
#     gca().clabel(CS, CS.levels, inline=true, fontsize=10)
#     xlabel("Latitude [°N]")
#     ylabel("Depth [km]")
#     gca().set_title(titlelabel)
#     gca().invert_yaxis()
#     gca().set_facecolor("black")
#     savefig(fname)
# end

# """
#     function sectionplotfancy
#     Plot of section (lat-depth) in ocean
# # Arguments
# - `field::Field`, 3d filed of values to be plotted
# - `lon`: longitude of section
# - `lims`: contour levels
# - `titlelabel`: optional title labeln
# """
# function sectionplotfancy(field::Field{<:Real}, lon, lims;titlelabel="section plot",fname="figure.png") 

#     Psection = section(field,lon)
#     cmap_seismic = get_cmap("seismic")
#     z = field.γ.depth/1000.0
    
#     #calc fignum - based on current number of figures
#     figure()
#     contourf(field.γ.lat, z, Psection', lims, cmap=cmap_seismic)
#     #fig, ax = plt.subplots()
#     CS = gca().contour(field.γ.lat, z, Psection', lims,colors="k")
#     gca().clabel(CS, CS.levels, inline=true, fontsize=10)
#     xlabel("Latitude [°N]")
#     ylabel("Depth [km]")
#     gca().set_title(titlelabel)
#     gca().invert_yaxis()
#     colorbar(orientation="horizontal")
#     savefig(fname)
# end

end
