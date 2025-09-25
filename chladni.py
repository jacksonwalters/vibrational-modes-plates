import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh
import argparse

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

# --- Main ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Chladni plate simulator")
    parser.add_argument("--bc", choices=["dirichlet", "neumann"], default="dirichlet",
                        help="Boundary condition type")
    parser.add_argument("--N", type=int, default=80, help="Grid points per axis")
    parser.add_argument("--modes", type=int, default=9, help="Number of modes to display")
    parser.add_argument("--skip", type=int, default=10, help="Number of lowest modes to skip")
    args = parser.parse_args()

    N = args.N
    L = 1.0
    num_modes = args.modes
    skip = args.skip
    bc = args.bc

    dx = L / (N - 1)
    L1D = laplacian_1d(N, dx, bc)
    n = L1D.shape[0]
    I = eye(n)
    L2D = kron(I, L1D) + kron(L1D, I)

    # Biharmonic operator
    B = L2D @ L2D

    # Eigenproblem
    vals, vecs = eigsh(B, k=num_modes+skip, sigma=0, which='LM')
    vals = vals[skip:]
    vecs = vecs[:, skip:]
    modes = [vecs[:, i].reshape((n, n)) for i in range(num_modes)]

    # Plotting
    X, Y = np.meshgrid(np.linspace(0, L, n), np.linspace(0, L, n))
    fig, axes = plt.subplots(3, 3, figsize=(9, 9))
    axes = axes.ravel()
    for i, ax in enumerate(axes):
        ax.contour(X, Y, modes[i], levels=[0], colors='k')
        ax.set_title(f"Mode {i+1}")
        ax.axis('off')

    plt.tight_layout()
    plt.show()
