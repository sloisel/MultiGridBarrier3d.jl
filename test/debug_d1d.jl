using LinearAlgebra
using SparseArrays

function chebyshev_nodes(k::Int)
    return [-cos(i * π / k) for i in 0:k]
end

function derivative_matrix_1d(nodes)
    k = length(nodes) - 1
    D = zeros(k+1, k+1)
    for i in 1:k+1
        xi = nodes[i]
        for j in 1:k+1
            if i == j
                val = 0.0
                for m in 1:k+1
                    if m != i
                        val += 1.0 / (xi - nodes[m])
                    end
                end
                D[i, j] = val
            else
                num = 1.0
                for m in 1:k+1
                    if m != j && m != i
                        num *= (nodes[i] - nodes[m])
                    end
                end
                den = 1.0
                for m in 1:k+1
                    if m != j
                        den *= (nodes[j] - nodes[m])
                    end
                end
                D[i, j] = num / den
            end
        end
    end
    return D
end

k = 3
nodes = chebyshev_nodes(k)
D = derivative_matrix_1d(nodes)

println("Nodes: $nodes")
println("D matrix:")
display(D)

# Test on linear function x
u = nodes
du = D * u
println("D * x (should be 1): $du")
println("Error: $(norm(du .- 1.0))")

# Test on quadratic x^2
u2 = nodes.^2
du2 = D * u2
println("D * x^2 (should be 2x): $du2")
println("Expected: $(2 .* nodes)")
println("Error: $(norm(du2 .- 2 .* nodes))")


# Test Tensor Product
I_1d = I(k+1)
D_xi = kron(I_1d, kron(I_1d, D))
D_eta = kron(I_1d, kron(D, I_1d))
D_zeta = kron(D, kron(I_1d, I_1d))

# Grid
nx = ny = nz = k+1
x = zeros(nx, ny, nz)
y = zeros(nx, ny, nz)
z = zeros(nx, ny, nz)

for iz in 1:nz, iy in 1:ny, ix in 1:nx
    x[ix, iy, iz] = nodes[ix]
    y[ix, iy, iz] = nodes[iy]
    z[ix, iy, iz] = nodes[iz]
end

x_flat = vec(x)
y_flat = vec(y)
z_flat = vec(z)

# Test D_eta on y
dy_deta = D_eta * y_flat
println("D_eta * y (should be 1): $(dy_deta[1:10])...")
println("Error: $(norm(dy_deta .- 1.0))")

# Test D_zeta on z
dz_dzeta = D_zeta * z_flat
println("D_zeta * z (should be 1): $(dz_dzeta[1:10])...")
println("Error: $(norm(dz_dzeta .- 1.0))")

# Test D_eta on x (should be 0)
dx_deta = D_eta * x_flat
println("D_eta * x (should be 0): $(norm(dx_deta))")
