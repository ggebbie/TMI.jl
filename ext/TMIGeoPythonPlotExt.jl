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
    function sectionplot
    Plot of section (lat-depth) in ocean
# Arguments
- `field::Field`, 3d filed of values to be plotted
- `lon`: longitude of section
- `lims`: contour levels
- `titlelabel`: optional title labeln
"""
function sectionplot(field::Field{<:Real}, lon, lims;titlelabel="section plot",fname="figure.png") 

    Psection = section(field,lon)

    GeoPythonPlot.sectionplot(Psection::Matrix, lon, field.γ.lat, field.γ.depth, lims;titlelabel="section plot",fname="figure.png",units=field.units) 
end

end
