```@meta
CurrentModule = MultiGridBarrier3d
```

```@eval
using Markdown
using Pkg
using MultiGridBarrier3d
v = string(pkgversion(MultiGridBarrier3d))
md"# MultiGridBarrier3d.jl $v"
```

A Julia package for solving 3D PDEs using the Spectral Barrier Method with hexahedral finite elements.

## Demo

Here is a simple example of solving a 3D problem and plotting the result.

```@example 1
using MultiGridBarrier3d

# Solve a 1-Laplacian problem on [-1,1]^3
sol = fem3d_solve(L=2, verbose=false)

# Plot the solution using PyVista
fig = plot(sol)
savefig(fig, "fem3d_demo.png"); nothing # hide
```

![](fem3d_demo.png)

## Installation

```julia
using Pkg
Pkg.add("MultiGridBarrier3d")
```
