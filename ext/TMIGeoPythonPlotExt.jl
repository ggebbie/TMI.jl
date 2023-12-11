module TMIGeoPythonPlotExt

using TMI
using GeoPythonPlot

#TMI.planviewplot()
"""
    function planviewplot: from NobleGasRelic

    formerly planviewplotcartopy
"""
function planviewplot(c::Field{<:Real}, depth, lims;titlelabel="planview plot",fname="fname.png",cenlon=-160.0) 

    cplan = planview(c,depth)

    # replace Field with Array
    GeoPythonPlot.planviewplot(cplan, depth, lims; titlelabel=titlelabel,fname=fname,cenlon=cenlon) 

end

#TMI.sectionplot()

end
