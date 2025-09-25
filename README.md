# vibrational-modes-plates

# Chladni Mode Simulation

This project simulates **Chladni figures** (the vibrational modes of plates and membranes) using finite differences in Python.

---

## Instructions

Clone the repository:

```sh
git clone https://github.com/jacksonwalters/vibrational-modes-plates.git
```

Create and activate a Python virtual environment called `math`:

```sh
python3 -m venv ~/.virtualenvs/math
source ~/.virtualenvs/math/bin/activate
```

Then simply run the script to produce plots:

```sh
python3 chlandi.py
```

## 1. Physical background

### Membrane (drum model)

A stretched membrane satisfies

\[
\rho\,\partial_{t}^{2}u = T\,\nabla^{2}u
\]

where:
- \(u(x,y,t)\) is displacement,
- \(\rho\) is surface density,
- \(T\) is tension.

Separation of variables \(u(x,y,t) = \phi(x,y)e^{i\omega t}\) gives

\[
\nabla^{2}\phi + \frac{\omega^{2}}{c^{2}}\,\phi = 0, \quad c^{2} = \tfrac{T}{\rho}.
\]

Equivalently,

\[
\nabla^{2}\phi = -k^{2}\phi, \qquad k^{2} = \tfrac{\omega^{2}}{c^{2}}.
\]

---

### Elastic plate (real Chladni plate)

A thin elastic plate obeys the Kirchhoff–Love equation:

\[
\rho h\,\partial_{t}^{2}w + D\,\nabla^{4}w = 0,
\]

with
- \(h\) plate thickness,
- \(D = \tfrac{E h^{3}}{12(1-\nu^{2})}\) bending stiffness.

Modes satisfy the **biharmonic eigenproblem**:

\[
\nabla^{4}\phi = \mu\,\phi, \quad \mu = \frac{\rho h}{D}\,\omega^{2}.
\]

---

## 2. Discrete Laplacian

On a 1-D grid with spacing \(\Delta x\):

\[
L_{1D} = \frac{1}{\Delta x^{2}}
\begin{bmatrix}
-2 & 1 & 0 & \cdots \\
1 & -2 & 1 & \ddots \\
0 & 1 & -2 & \ddots \\
\vdots & \ddots & \ddots & \ddots
\end{bmatrix}.
\]

On a 2-D square grid:

\[
L_{2D} = I \otimes L_{1D} + L_{1D} \otimes I.
\]

This matrix approximates \(\nabla^{2}\phi\).

- Eigenpairs satisfy \(\nabla^{2}\phi = -k^{2}\phi\).
- Numerically: if \(L_{2D} v = \lambda v\), then \(k^{2} \approx -\lambda\).

---

## 3. Boundary conditions

- **Dirichlet (clamped edge)**: enforce \(\phi = 0\) on the boundary (remove or fix boundary rows/cols).
- **Neumann (free edge)**: enforce \(\partial_n \phi = 0\) (zero normal derivative at the boundary).
- **Plate BCs**: more complex (need both displacement and slope/moment conditions).

---

## 4. Eigenvalues → frequencies

- **Membrane**: \(-L_{2D} v = k^{2} v \implies \omega = c k\).
- **Plate**: \((L_{2D})^{2} v = \mu v \implies \omega = \sqrt{\tfrac{D}{\rho h}\mu}\).

---

## 5. Chladni patterns

Nodal lines are given by the zero set of \(\phi(x,y)\).  
In experiments, sand accumulates along these nodal lines.  
Numerically: plot contour lines at level = 0.

---

## 6. Implementation notes

- Use `scipy.sparse` and `eigsh` for larger grids.
- For membrane: solve eigenpairs of `-L2D`.
- For plate: approximate \(\nabla^{4}\) with `(L2D)^2` (basic clamped-edge model).
- Plot zero contours for Chladni figures.

