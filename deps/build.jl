import Pkg; Pkg.add("Pkg")

using Pkg

if lowercase(get(ENV, "CI", "false")) == "true"    

        ENV["PYTHON"] = ""
        Pkg.build("PyCall")
      
end
