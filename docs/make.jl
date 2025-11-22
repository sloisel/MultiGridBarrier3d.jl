using Pkg
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using MultiGridBarrier3d
using PyPlot

makedocs(
    sitename = "MultiGridBarrier3d.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [MultiGridBarrier3d],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ],
    # remotes = nothing,
    checkdocs = :exports
)

deploydocs(
    repo = "github.com/sloisel/MultiGridBarrier3d.jl.git",
    devbranch = "main",
)
