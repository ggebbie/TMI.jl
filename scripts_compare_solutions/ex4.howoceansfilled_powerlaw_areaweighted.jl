#=
% Example: Find the ocean volume that has originated from each    %
%            surface box.                                           %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of each box %
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.
%
% See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011. 
=#

import Pkg; Pkg.activate("./scripts_compare_solutions")
Pkg.add("PythonPlot")
Pkg.instantiate()

# using Revise
using TMI
using Test
using GeoPythonPlot
using PythonPlot

TMI_versions = TMI.versionlist()[1:end-1]

fig, ax = subplots(2, 5, figsize = (10, 7.5), sharey = true, sharex = true)
ax_flat = ax.flatten()
for (i, TMIversion) in enumerate(TMI_versions)
    ax_flat[i-1].set_title(TMIversion) #i-1 since its a Python object
    ax_flat[i-1].set_xlabel("cumulative area")
    ax_flat[i-1].set_ylabel("log₁₀(weight)")

    try
        A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
        areas = cellarea(γ)
        volume = volumefilled(TMIversion,Alu,γ)
        sort_idx = sortperm(volume.tracer[:], rev = true)
        sorted_weights = volume.tracer[:][sort_idx]
        sorted_areas = areas.tracer[:][sort_idx]
        cumulative_areas = cumsum(sorted_areas)
        ax_flat[i-1].scatter(cumulative_areas, sorted_weights)
    catch
        print("Can't grab TMIVersion: ", TMIversion)
        print("Probably an issue with wget")
    end
end
fig.tight_layout()
fig