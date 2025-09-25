import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# Parameters
N = 50   # grid resolution
L = 1.0  # plate length (square domain)

# Discretization
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
dx = x[1] - x[0]
dy = y[1] - y[0]
X, Y = np.meshgrid(x, y, indexing="ij")

# Laplacian in 2D with Dirichlet boundary conditions
def laplacian_2d(N, dx):
    """Construct 2D Laplacian with Dirichlet BC on N x N grid."""
    I = np.eye(N)
    D = np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1),1) + np.diag(np.ones(N-1),-1)
    L1D = D / dx**2
    # Kronecker sum
    L2D = np.kron(I, L1D) + np.kron(L1D, I)
    return L2D

# Build Laplacian
L2D = laplacian_2d(N, dx)

# Solve eigenvalue problem: L u = λ u
# We want the lowest few modes
num_modes = 9
vals, vecs = eigh(L2D)

# Pick the lowest positive eigenmodes
modes = []
for i in range(num_modes):
    u = vecs[:, i].reshape(N, N)
    modes.append(u)

# Plot results
fig, axes = plt.subplots(3, 3, figsize=(9, 9))
for ax, mode, val in zip(axes.ravel(), modes, vals[:num_modes]):
    im = ax.contourf(X, Y, mode, levels=20, cmap="RdBu")
    ax.set_title(f"λ = {val:.2f}")
    ax.axis("off")
plt.tight_layout()
plt.show()
