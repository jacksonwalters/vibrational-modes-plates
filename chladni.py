import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, eye
from scipy.sparse.linalg import eigsh

# Parameters
N = 80            # grid points per axis (higher = finer patterns)
L = 1.0
num_modes = 9     # number of modes to display
skip = 10         # number of lowest modes to skip

# --- 1D Laplacian (Dirichlet BCs) ---
def laplacian_1d(N, dx):
    n = N - 2
    main = -2 * np.ones(n)
    off = np.ones(n - 1)
    return diags([off, main, off], [-1, 0, 1]) / dx**2

# --- Build 2D Laplacian ---
n = N - 2
dx = L / (N - 1)
L1D = laplacian_1d(N, dx)
I = eye(n)
L2D = kron(I, L1D) + kron(L1D, I)

# --- Biharmonic operator (plate) ---
B = L2D @ L2D

# --- Solve eigenproblem ---
# Use shift-invert to get small eigenvalues if needed, or eigsh normally
vals, vecs = eigsh(B, k=num_modes+skip, sigma=0, which='LM')
vals = vals[skip:]
vecs = vecs[:, skip:]
modes = [vecs[:, i].reshape((n, n)) for i in range(num_modes)]

# --- Plot ---
X, Y = np.meshgrid(np.linspace(0, L, n), np.linspace(0, L, n))
fig, axes = plt.subplots(3, 3, figsize=(9, 9))
axes = axes.ravel()
for i, ax in enumerate(axes):
    ax.contour(X, Y, modes[i], levels=[0], colors='k')
    ax.set_title(f"Mode {i + 1}")
    ax.axis('off')

plt.tight_layout()
plt.show()
