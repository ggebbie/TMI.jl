# TMI

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/TMI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/TMI.jl/dev)
[![Build Status](https://github.com/ggebbie/TMI.jl/workflows/CI/badge.svg)](https://github.com/ggebbie/TMI.jl/actions)
[![Coverage](https://codecov.io/gh/ggebbie/TMI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/TMI.jl)

* Total Matrix Intercomparison codes for Julia\
Started by G Jake Gebbie, WHOI, ggebbie@whoi.edu 

* See the function list in the documentation linked through the badge above

* The MATLAB version of the code is in maintenance mode and is available at https://github.com/ggebbie/TMI 

* After setting up the environment (instructions below), check that all tests pass via the following shell command in the repository base directory:
`julia --project=@. test/runtests.jl`

# Requirements

The built-in tests are automatically checked with Julia 1.8. 

# Setting up project environment

* from the Julia REPL

`;`\
`git clone https://github.com/ggebbie/TMI.jl # only do this the first time on each machine`\
`cd TMI.jl`\
`]`\
`activate .`\
`instantiate # only do this the first time on each machine`\
To verify you are in the project environment, `]` should return `(TMI) pkg>`\
Type backspace to return to command mode.

* from Emacs editor (one possible method)

Install julia-mode, julia-repl or julia-snail, and magit \
Skip the next 5 steps if you have already cloned the repository \
`M-x magit-clone` \
Select `u` to clone from url\
Enter ` https://github.com/ggebbie/TMI.jl` as url to clone \
Select `y` in response to `remote.pushDefault' to "origin"?` \
Clone to your favorite location and rename project if necessary \
Go to any directory in the project: `C-x C-f TMI.jl`\
Then activate the project and initialize a julia session: `C-c C-a`

* Using an editor like Atom/Juno or Visual Studio Code, activate the environment on one of the frame panels. The default environment is @v1.x and should be changed.

# Scripts/Examples

Examples for what TMI can do are found in the `scripts` directory.

See examples in `scripts/ex*`, where `ex1.trackpathways.jl` gives Example 1 of tracking water-mass pathways, for example. Scripts beginning with "ex" are tested and can serve as good templates for new work. Other scripts are a work in progress.

Scripts can be run non-interactively like this:\
`cd TMI.jl`\
`julia --project=scripts scripts/ex1.trackpathways.jl`

To show graphical results, TMI.jl uses `GGplot.jl` for plotting routines. In particular, matplotlib, cartopy and cmocean packages are handled by `GGplot` so that they are not dependencies in `TMI.jl`. Thus the `scripts` directory has its own environment distinct from the TMI project. If you are working interactively, try the following commands to activate the scripts environment: 

`cd("scripts")`\
`import Pkg; Pkg.activate(".")`

If you use the examples and `GGplot`, you will need a python environment in Julia. Direct the python environment to an existing system-wide version of python with these already installed:
`ENV["PYTHON"]="python/directory/on/your/machine"`

GGplot will use a package-specific python environment built from scratch using the `CondaPkg.jl` package. Check out the `GGplot/deps/build.jl` file to see how this Python environment is set up. In particular, it executes:
`ENV["PYTHON"]=""` 

Rather than using the `PyCall.jl` package, `GGplot.jl` uses `PythonCall.jl` in order to minimize errors that occur due to incompatible Python environments. 
`import Pkg; Pkg.add("PythonCall")`\
`import Pkg; Pkg.build("PythonCall") # probably not necessary`

In particular, GGplot uses the matplotlib, cartopy, and cmocean packages. The channel is automatically assumed to be conda-forge. \
`using CondaPkg`\
`] conda add matplotlib # from the package manager`\
`] Conda add cartopy`\
`] Conda add cmocean`

This version of cartopy (v0.20.0+) can download coastlines from Amazon Web Services.

# Data files

The Julia code is designed to download input files from Google Drive and to place them in the `data` directory. If that doesn't work, extract data from Google Drive using your favorite method or download manually at: https://drive.google.com/drive/folders/1nFDSINbst84pK68aWwRGBVfYZkfN1mUR?usp=sharing . Available TMI versions include:

`modern_90x45x33_GH10_GH12` : TMI version with 4x4 degree horizontal
                  resolution and 33 levels  (G & H 2010), \
				  Includes the input data from the WGHC (Gouretski & Koltermann 2005) 
 
`modern_180x90x33_GH11_GH12` : TMI version with 2x2 degree horizontal
                  resolution and 33 levels  (G & H 2011), \
				  Includes the input data from the WGHC (Gouretski & Koltermann 2005) 

`modern_90x45x33_unpub12` : TMI version with 4x4 degree horizontal
                  resolution and 33 levels  (unpublished 2012), \
				  Includes a steady-state climatology of global tracers

`modern_90x45x33_G14` : TMI version with 4x4 degree horizontal
                  resolution and 33 levels  (Gebbie 2014), \
				  Doesn't rely upon a bottom spreading parameterization and solves for mixed-layer depth

`modern_90x45x33_G14_v2` : TMI version with 4x4 degree horizontal
                  resolution and 33 levels  (Gebbie 2014), \
				  Doesn't rely upon a bottom spreading parameterization and solves for mixed-layer depth\
				  Includes optimization information
				  
`LGM_90x45x33_G14` : Last Glacial Maximum version with 4x4 degree horizontal
                  resolution and 33 levels  (Gebbie 2014)
				  
`LGM_90x45x33_G14A` : Alternate solution, Last Glacial Maximum version with 4x4 degree horizontal
                  resolution and 33 levels  (Gebbie 2014)
				  
`LGM_90x45x33_GPLS1`: Solution #1 (Gebbie, Peterson, Lisiecki, and Spero, 2015), Last Glacial Maximum version with 4x4 degree horizontal resolution and 33 levels 
				  
`LGM_90x45x33_GPLS2`: Solution #2 (Gebbie, Peterson, Lisiecki, and Spero, 2015), Last Glacial Maximum version with 4x4 degree horizontal resolution and 33 levels 
				  
`LGM_90x45x33_OG18`: Last Glacial Maximum version with 4x4 degree horizontal resolution and 33 levels (Oppo, Gebbie et al. 2018)

`nordic_201x115x46_B23`: Nordic TMI regional version, 1/2 degree longitude by 1/4 degree latitude and 46 levels (Brakstad et al., 2013) 

# Functions

Available functions are listed in the documentation at https://ggebbie.github.io/TMI.jl/dev/ .


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

# References 

Gebbie, G., and P. Huybers:  "Total matrix intercomparison: A method for resolving the geometry of water mass pathways", J. Phys. Oceanogr., 40(8), doi:10.1175/2010JPO4272.1, 1710-1728, 2010. 

Gebbie, G., and P. Huybers. "How is the ocean filled?", Geophys. Res. Lett., 38, L06604, doi:10.1029/2011GL046769, 2011 

Gebbie, G., and P. Huybers, "The mean age of ocean waters inferred from radiocarbon observations", 2012, JPO.

Gebbie, G., "How much did Glacial North Atlantic Water shoal?", 2014, Paleoceanography.

# How this Julia package was started

This package was generated using PkgTemplates.jl. 

Steps: 
1. Use PkgTemplates to make git repo.\
	2. new empty repository on GitHub.\
    3. Then push an existing repository from the command line:
    `git remote add origin git@github.com:ggebbie/TMI.jl.git`\
    `git branch -M main`\
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
        GitHubActions(; x86=false),
        Codecov(),
        Documenter{GitHubActions}(),
        Develop(),
    ],
             )`

`t("TMI.jl")`
