using Revise, TMI
using Plots
#= 
 Example 1: Track the pathways of a user-defined water mass.         
 Steps: (a) define the water mass 1). by a rectangular surface patch 
            dyed with passive tracer concentration of 1,              
            or 2. load a pre-defined surface patch in d_all.mat.     
        (b) propagate the dye with the matrix A, with the result      
            being the fraction of water originating from the         
            surface region.                                          
 See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).       
=#

url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
inputdir = "../input"

A, Alu, γ = configTMI(url,inputdir)

#- define the surface patch by the bounding latitude and longitude.
lat_lo = 50; # 50 N, for example.
lat_hi = 60;
 
lon_lo = -50;
lon_hi = 0;

# ternary operator to handle longitudinal wraparound
lon_lo ≤ 0 ? lon_lo += 360 : nothing
lon_hi ≤ 0 ? lon_hi += 360 : nothing

# define the surface boundary condition
nfield = size(A,1)
d = zeros(Int64,nfield) # preallocate 

# Extract i, j, and k indices (must be a better way to do this)
i = zeros(Int,nfield)
j = zeros(Int,nfield)
k = zeros(Int,nfield)
d = zeros(Real,nfield)

[i[n] = γ.coords[n][1] for n ∈ 1:nfield]
[j[n] = γ.coords[n][2] for n ∈ 1:nfield]
[k[n] = γ.coords[n][3] for n ∈ 1:nfield]
    
[d[n]=1 for n ∈ 1:nfield if k[n]==1 && lat_lo ≤ γ.lat[j[n]] ≤ lat_hi && lon_lo ≤ γ.lon[i[n]] ≤ lon_hi ]

# do matrix inversion to get quantity of dyed water throughout ocean:
c = Alu\d

# after doing calculations with vectors, translate to a 3D geometric field
cfld = vec2fld(c,γ.coords)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330;
isec = findall(==(lon_section),γ.lon)
csection = dropdims(cfld[isec,:,:],dims=1)

# make a plot of dye in the ocean
Plots.contourf(γ.lat,-γ.depth[33:-1:1],csection[:,33:-1:1]') # a sample plot at 22 W.
