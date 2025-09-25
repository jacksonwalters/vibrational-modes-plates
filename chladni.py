import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh

# --- Parameters ---
N = 80
L = 1.0
num_modes = 9
skip = 10
bc = "neumann"  # "dirichlet" or "neumann"

# --- 1D Laplacian ---
def laplacian_1d(N, dx, bc="dirichlet"):
    if bc == "dirichlet":
        n = N - 2
        main = -2 * np.ones(n)
        off = np.ones(n - 1)
        return diags([off, main, off], [-1, 0, 1]) / dx**2
    elif bc == "neumann":
        n = N
        main = -2 * np.ones(n)
        off = np.ones(n - 1)
        L = diags([off, main, off], [-1, 0, 1]) / dx**2
        L = L.tolil()
        L[0,0] = -1/dx**2; L[0,1] = 1/dx**2
        L[-1,-1] = -1/dx**2; L[-1,-2] = 1/dx**2
        return csr_matrix(L)
    else:
        raise ValueError("bc must be 'dirichlet' or 'neumann'")

dx = L / (N - 1)
L1D = laplacian_1d(N, dx, bc)

n = L1D.shape[0]
I = eye(n)
L2D = kron(I, L1D) + kron(L1D, I)

# --- Biharmonic operator ---
B = L2D @ L2D

# --- Solve eigenproblem ---
vals, vecs = eigsh(B, k=num_modes+skip, sigma=0, which='LM')
vals = vals[skip:]
vecs = vecs[:, skip:]
modes = [vecs[:, i].reshape((n, n)) for i in range(num_modes)]

# --- Plotting ---
X, Y = np.meshgrid(np.linspace(0, L, n), np.linspace(0, L, n))
fig, axes = plt.subplots(3, 3, figsize=(9, 9))
axes = axes.ravel()
for i, ax in enumerate(axes):
    ax.contour(X, Y, modes[i], levels=[0], colors='k')
    ax.set_title(f"Mode {i+1}")
    ax.axis('off')

plt.tight_layout()
plt.show()
