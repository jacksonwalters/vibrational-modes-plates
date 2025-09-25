# vibrational-modes-plates

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

# Chladni Mode Simulation

This project simulates **Chladni figures** (the vibrational modes of plates) using finite differences in Python.

---

## 1. Physical background

### Elastic plate (Chladni plate)

A thin elastic plate obeys the Kirchhoff–Love equation:

$$
\rho h\,\partial_{t}^{2}w + D\,\nabla^{4}w = 0,
$$

with  
- $h$ = plate thickness  
- $D = \frac{E h^{3}}{12(1-\nu^{2})}$ = bending stiffness

Modes satisfy the **biharmonic eigenproblem**:

$$
\nabla^{4}\phi = \mu\,\phi, \quad \mu = \frac{\rho h}{D}\,\omega^{2}.
$$

---

## 2. Discrete Laplacian

On a 1-D grid with spacing $\Delta x$:

$$
L_{1D} = \frac{1}{\Delta x^{2}}
\begin{bmatrix}
-2 & 1 & 0 & \cdots \\
1 & -2 & 1 & \ddots \\
0 & 1 & -2 & \ddots \\
\vdots & \ddots & \ddots & \ddots
\end{bmatrix}.
$$

On a 2-D square grid:

$$
L_{2D} = I \otimes L_{1D} + L_{1D} \otimes I.
$$

This matrix approximates $\nabla^{2}\phi$. For a plate, the biharmonic operator is $L_{2D}^2$.

---

## 3. Boundary conditions

- **Dirichlet (clamped edge)**: enforce $\phi = 0$ on the boundary.  
- **Neumann (free edge)**: enforce $\partial_n \phi = 0$ (zero normal derivative).  

---

## 4. Eigenvalues → frequencies

- **Plate**: $(L_{2D})^2 v = \mu v \implies \omega = \sqrt{\frac{D}{\rho h}\mu}$.

---

## 5. Chladni patterns

Nodal lines are given by the zero set of $\phi(x,y)$.  
Numerically: plot contour lines at level = 0 to simulate where sand would accumulate.

---

## 6. Implementation notes

- Use `scipy.sparse` and `eigsh` for efficiency.  
- For realistic clamped edges, enforce Dirichlet BCs (interior points only).  
- Plot zero contours (`plt.contour(..., levels=[0])`) to visualize Chladni patterns.  
- Increase grid resolution for finer patterns.  
- Skip the first few lowest modes to see higher-order “flower” and circular patterns.


