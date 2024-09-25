using Documenter
using PML

if haskey(ENV, "GITHUB_ACTIONS")
  ENV["JULIA_DEBUG"] = "Documenter"
end

deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type="pending", repo="github.com/jl-pkgs/PML.jl.git")
using Literate
using Plots # to not capture precompilation output

# generate examples
EXAMPLE = joinpath(@__DIR__, "..", "examples", "CalibOneSite.jl")
OUTPUT = joinpath(@__DIR__, "src/generated")

preprocess(str) = replace(str, "x = 123" => "y = 321"; count=1)

# Literate.notebook(EXAMPLE, OUTPUT; preprocess)
Literate.markdown(EXAMPLE, OUTPUT; preprocess)
# Literate.script(EXAMPLE, OUTPUT; preprocess)

# # generate the example notebook for the documentation, keep in sync with outputformats.md
# Literate.markdown(joinpath(@__DIR__, "src/outputformats.jl"), OUTPUT; credit=false, name="name")
# Literate.notebook(joinpath(@__DIR__, "src/outputformats.jl"), OUTPUT; name="notebook")
# Literate.script(joinpath(@__DIR__, "src/outputformats.jl"), OUTPUT; credit=false)

# Replace the link in outputformats.md
# since that page is not "literated"
if haskey(ENV, "GITHUB_ACTIONS")
  folder = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
    Documenter.deploy_folder(
      deployconfig;
      repo="github.com/jl-pkgs/PML.jl.git",
      devbranch="master",
      push_preview=true,
      devurl="dev",
    ).subfolder
  end
  url = "https://nbviewer.jupyter.org/github/jl-pkgs/PML.jl/blob/gh-pages/$(folder)/"
  # str = read(joinpath(@__DIR__, "src/outputformats.md"), String)
  # str = replace(str, "[notebook.ipynb](generated/notebook.ipynb)." => "[notebook.ipynb]($(url)generated/notebook.ipynb).")
  # write(joinpath(@__DIR__, "src/outputformats.md"), str)
end

# # Generate changelog
# using Changelog
# clog = joinpath(@__DIR__, "src", "changelog.md")
# Changelog.generate(
#   Changelog.Documenter(),
#   joinpath(@__DIR__, "..", "CHANGELOG.md"),
#   clog;
#   repo="jl-pkgs/PML.jl",
# )
# write(clog, replace(read(clog, String), r"^# PML.jl changelog"m => "# **9.** Changelog"))


makedocs(
  format=Documenter.HTML(
    # assets=["assets/custom.css", "assets/favicon.ico"],
    prettyurls=true, # haskey(ENV, "GITHUB_ACTIONS"),
    size_threshold=1024*1024,
    # canonical="https://jl-pkgs.github.io/PML.jl",
  ),
  modules=[PML],
  sitename="PML.jl",
  pages=Any[
    "index.md",
    # "fileformat.md",
    # "pipeline.md",
    # "outputformats.md",
    # "customprocessing.md",
    # "documenter.md",
    # "tips.md",
    "generated/CalibOneSite.md",
    # "changelog.md",
  ],
  warnonly=true,
)

deploydocs(
  repo="github.com/jl-pkgs/PML.jl.git",
  push_preview=true,
  versions=["v2" => "v^", "v#.#", "dev" => "dev"],
  deploy_config=deployconfig,
)
