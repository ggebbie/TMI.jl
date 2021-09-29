# TMI

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/TMI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/TMI.jl/dev)
[![Build Status](https://github.com/ggebbie/TMI.jl/workflows/CI/badge.svg)](https://github.com/ggebbie/TMI.jl/actions)
[![Coverage](https://codecov.io/gh/ggebbie/TMI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ggebbie/TMI.jl)

* Total Matrix Intercomparison codes for Julia\
G Jake Gebbie, WHOI, ggebbie@whoi.edu 

# References 

Gebbie, G., and P. Huybers:  "Total matrix intercomparison: A method for resolving the geometry of water mass pathways", J. Phys. Oceanogr., 40(8), doi:10.1175/2010JPO4272.1, 1710-1728, 2010. 

Gebbie, G., and P. Huybers. "How is the ocean filled?", Geophys. Res. Lett., 38, L06604, doi:10.1029/2011GL046769, 2011 

Gebbie, G., and P. Huybers, "The mean age of ocean waters inferred from radiocarbon observations", 2012, JPO.

Gebbie, G., "How much did Glacial North Atlantic Water shoal?", 2014, Paleoceanography.

# Requirements

TMI.jl uses matplotlib and cartopy. Direct the python environment to an existing system-wide version of python with these already installed:
`ENV["PYTHON"]="python/directory/on/your/machine"`

Or use a Julia-specific python environment built from scratch following these directions from the Julia REPL:
`ENV["PYTHON"]=""` \
`import Pkg; Pkg.add("PyCall")`\
`import Pkg; Pkg.build("PyCall")`\
Restart the REPL.\
`import Pkg; Pkg.add("Conda")`\
`import Conda`\
`Conda.add("matplotlib",channel="conda-forge")`\
`Conda.add("shapely",channel="conda-forge")`\
`Conda.add("cartopy",channel="conda-forge")`\

This should set up cartopy v. 0.20.0 which can download coastlines from Amazon Web Services.

# Setting up project environment

* from Emacs editor (one possible method)

Install julia-mode, julia-repl, and magit \
Skip the next 5 steps if you have already cloned the repository \
`M-x magit-clone` \
Select `u` to clone from url\
Enter ` https://github.com/ggebbie/TMI.jl` as url to clone \
Select `y` in response to `remote.pushDefault' to "origin"?` \
Clone to your favorite location and rename project if necessary \
Go to any directory in the project: `C-x C-f TMI.jl`\
Then activate the project and initialize a julia session: `C-c C-a`

* from the Julia REPL

`;`\
`git clone https://github.com/ggebbie/TMI.jl # only do this the first time on each machine`\
`cd TMI.jl`\
`]`\
`activate .`\
`instantiate # only do this the first time on each machine`\
To verify you are in the project environment, `]` should return `(TMI) pkg>`\
Type backspace to return to command mode.

* Using an editor like Atom/Juno or Visual Studio Code, activate the environment on one of the frame panels. The default environment is @v1.x and should be changed.

# Running a script (not interactively)

An example:\
`cd TMI.jl`\
`julia --project=@. scripts/ex1.trackpathways.jl`

# Data files

The Julia code is designed to download input files from Google Drive and to place them in the `data` directory. If that doesn't work, extract data from Google Drive using your favorite method, such as the script `readTMIfromGoogleDrive.sh` at a bash shell prompt. 

`/bin/sh readTMIfromGoogleDrive.sh`

Or download manually at: https://drive.google.com/file/d/1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7/view?usp=sharing .

TMI_4deg_2010.mat : TMI version with 4x4 degree horizontal
                  resolution and 33 levels  (G & H 2010), \
				  Includes TMI climatology of ocean properties 

# Functions

Available functions are listed in the documentation at https://ggebbie.github.io/TMI.jl/dev/ .

# Scripts

See examples in `scripts/ex*`, where `ex1.trackpathways.jl` gives Example 1 of tracking water-mass pathways, for example.

# How this Julia package was started

This package was generated using PkgTemplates.jl. 

Steps: 
1. Use PkgTemplates to make git repo.\
	2. new empty repository on GitHub.\
    3. Then push an existing repository from the command line:
    `git remote add origin git@github.com:ggebbie/TMI.jl.git`
    `git branch -M main`
    `git push -u origin main`

4. Run the following Julia code

`using Revise, PkgTemplates`

`t = Template(; 
    user="ggebbie",
    dir="~/projects",
    authors="G Jake Gebbie",
    julia=v"1.6",
    plugins=[
        License(; name="MIT"),
        Git(; manifest=true, ssh=true),
        GitHubActions(; x86=true),
        Codecov(),
        Documenter{GitHubActions}(),
        Develop(),
    ],
             )`

`t("TMI.jl")`

# MATLAB version of code

MATLAB codes, 2009-2021, see also https://github.com/ggebbie/TMI .

History:\
Version 1, 7 May 2009.\
Version 2, 6 Aug 2010.\
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
