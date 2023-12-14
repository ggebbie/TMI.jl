module TMIGeoPythonPlotExt

using TMI
using GeoPythonPlot
import GeoPythonPlot: planviewplot, sectionplot

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

"""
    function planviewplot
    Plot of plan view (lon-lat) in ocean
# Arguments
- `field::BoundaryCondition`, 3d field of values to be plotted
- `lims`: contour levels
# Optional Arguments
- `titlelabel`: optional title label
- `fname`: file name for saving plot
- `cenlon`: central longitude
"""
function planviewplot(b::BoundaryCondition, lims; titlelabel="boundary condition plot",fname="figure.png", cenlon= -160.0)

    lon = b.i
    lat = b.j
    depth = b.k
    
    # replace Field with Array
    GeoPythonPlot.planviewplot(b.tracer, lon, lat, depth, lims; titlelabel=titlelabel,fname=fname,cenlon=cenlon,units=b.units) 

end

#     # is the boundary condition oriented correctly?
#     if b.dim != 3
#         error("boundary condition not horizontal")
#     end
    
#     cplan = b.tracer
    
#     cmap_seismic = get_cmap("seismic")
    
#     #calc fignum - based on current number of figures
#     figure()
#     contourf(γ.lon,γ.lat, cplan', lims, cmap=cmap_seismic)
#     #fig, ax = plt.subplots()
#     CS = gca().contour(γ.lon,γ.lat, cplan', lims, cmap=cmap_seismic)
#     gca().clabel(CS, CS.levels, inline=true, fontsize=10)
#     ylabel("Latitude [°N]")
#     xlabel("Longitude [°E]")
#     gca().set_title(titlelabel)
#     colorbar(orientation="vertical")
    
# end

"""
    function sectionplot
    Plot of section (lat-depth) in ocean
# Arguments
- `field::Field`, 3d filed of values to be plotted
- `lon`: longitude of section
- `lims`: contour levels
- `titlelabel`: optional title labeln
"""
function sectionplot(field::Field{<:Real}, lon, lims;titlelabel="section plot",fname="figure.png",units=:none) 

    Psection = section(field,lon)

    GeoPythonPlot.sectionplot(Psection::Matrix, lon, field.γ.lat, field.γ.depth, lims;titlelabel="section plot",fname="figure.png",units=field.units) 
end

end
