using Pkg
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using MultiGridBarrier3d
using PyPlot

# Compute version dynamically
version = string(pkgversion(MultiGridBarrier3d))

makedocs(
    sitename = "MultiGridBarrier3d.jl $version",
    remotes = nothing,
    warnonly = [:missing_docs, :cross_references, :docs_block],
    modules = [MultiGridBarrier3d],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://sloisel.github.io/MultiGridBarrier3d.jl/stable/",
    ),
    checkdocs = :exports
)

deploydocs(
    repo = "github.com/sloisel/MultiGridBarrier3d.jl",
    devbranch = "main",
)
