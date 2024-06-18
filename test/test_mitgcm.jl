@testset "mitgcm" begin

    using GoogleDrive
    
    # construct Grid from MITgcm output
    url = "https://docs.google.com/uc?export=download&id=16-F3V-MUi3cIOqgJiJ4x4JD0sGlHGQ1y"
    fname = google_download(url,TMI.pkgdatadir())
    println(fname)
    γ =  Grid(fname, "maskC", "XC", "YC", "Z")

    # Add a more stringest test here, like a `checkgrid!`

    θ = readfield(fname,"THETA",γ)
end
