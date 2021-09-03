# TMI

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/TMI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/TMI.jl/dev)
[![Build Status](https://github.com/ggebbie/TMI.jl/workflows/CI/badge.svg)](https://github.com/ggebbie/TMI.jl/actions)
[![Coverage](https://codecov.io/gh/ggebbie/TMI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ggebbie/TMI.jl)

* Total Matrix Intercomparison MATLAB diagnostic routines 
Codes by Geoffrey (Jake) Gebbie, WHOI, ggebbie@whoi.edu, 2009-2021

* References 

Gebbie, G., and P. Huybers:  "Total matrix intercomparison: A method for resolving the geometry of water mass pathways", J. Phys. Oceanogr., 40(8), doi:10.1175/2010JPO4272.1, 1710-1728, 2010. 

Gebbie, G., and P. Huybers. "How is the ocean filled?", Geophys. Res. Lett., 38, L06604, doi:10.1029/2011GL046769, 2011 

Gebbie, G., and P. Huybers, "The mean age of ocean waters inferred from radiocarbon observations", 2012, JPO.

Gebbie, G., "How much did Glacial North Atlantic Water shoal?", 2014, Paleoceanography.

* History of MATLAB codes

Version 1, 07 May 2009.\
Version 2, 06 Aug 2010.\
Version 3, 21 Apr 2011 -- minor changes.\
Version 4, 13 July 2011, makes names consistent with papers.\
Version 5, 28 July 2011, add TMI transient tracer simulation model.\
Version 6, Nov 2012, bug fixes, use one LU decomp for both fwd and
                        adjoint, added global inversion example,
                        SynTraCE-21 workshop update \
Version 6.1, Jan 2013, added biogeochemical example, add
                       vector_to_field back into tarball.\
Version 6.2, July 2015, added sq.m function,
                        fixed d_all to properly divide Atlantic/Pacific and put White Sea into Arctic.\
Version 7, Sept. 2016, major improvements to transient run: 2 types of initial conditions and boundary conditions.\
Version 8, Jan. 2021, bug fixes, especially those found by Elaine McDonagh

# MAIN DIAGNOSTIC ROUTINES:

Directory `scripts`: \ 
`steadystate_diagnostics.m`  : examples of analysis for the TMI pathways matrix.\
`transient_driver.m` : run a TMI transient tracer simulation model.

# DATA FILES

Extract data from Google Drive using your favorite method. Try the script `read_TMI_from_google_drive` either in MATLAB or at a bash shell prompt. 

`/bin/sh read_TMI_from_google_drive.sh`

Or download manually at: https://drive.google.com/file/d/1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7/view?usp=sharing .

TMI_4deg_2010.mat : TMI version with 4x4 degree horizontal
                  resolution and 33 levels  (G & H 2010) \
				  Includes TMI climatology of ocean properties \

# FUNCTIONS 

Source code in  `src/TMI.jl`. \ 
See the list of `export`ed functions in the header of that file.
