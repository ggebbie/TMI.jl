""" 
    A short code to show how this package was made
    according to a template. This is a utility script that
    was run just once to create the package and will not
    be needed for this project again.

    Steps: 
    1. Use PkgTemplates to make git repo.
    2. new empty repository on GitHub.
    3. Then push an existing repository from the command line:
    `git remote add origin git@github.com:ggebbie/TMI.jl.git`
    `git branch -M main`
    `git push -u origin main`

"""
using Revise, PkgTemplates

t = Template(; 
    user="ggebbie",
    dir="~/projects",
    authors="G Jake Gebbie",
    julia=v"1.9",
    plugins=[
        License(; name="MIT"),
        Git(; manifest=true, ssh=true),
        GitHubActions(; x86=true),
        Codecov(),
        Documenter{GitHubActions}(),
        Develop(),
    ],
             )

t("TMI.jl")
