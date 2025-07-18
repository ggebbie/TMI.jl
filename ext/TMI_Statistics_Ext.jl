module TMI_Statistics_Ext

using TMI, Statistics

"""
    function Statistics.mean(c::Field)

    Take the volume-weighted mean of a `Field`
"""
function Statistics.mean(c::Field)
    vol = cellvolume(c.γ)
    return sum(vol.tracer[wet(c)].*c.tracer[wet(c)])/sum(vol.tracer[wet(c)])
end


end
