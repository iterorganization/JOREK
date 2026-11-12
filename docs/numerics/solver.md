---
title: "Solver and Preconditioner"
nav_order: 5
parent: "Numerics and Tools"
layout: default
render_with_liquid: false
---

# JOREK Solver
The implicit time integration scheme implemented by JOREK leads to a **large** and **sparse** linear system of the general form

$$A\,x = b.$$

The matrix $A$ is assembled in a distributed block-COO format (which can be converted to a distributed block-CSR format on demand) and is partitioned by **toroidal Fourier harmonics** and **finite-element nodes**. The solver infrastructure exploits this block structure at both the factorisation and the preconditioning level.

Two high-level solution strategies are supported:

| | **Direct solve** | **Iterative solve** |
|---|---|---|
| Input flag | `gmres = .false.` | `gmres = .true.` |
| Method | Sparse LU factorisation of $A$ | GMRES (or BiCGSTAB) with block preconditioner |
| Strengths | Robust; insensitive to conditioning | Lower memory; scales to large problems |
| Weaknesses | Memory and fill-in grow rapidly | Convergence depends on preconditioner |
| Typical use | Small-to-moderate problems | Large / production runs |

Both strategies support multiple third-party libraries as back-ends and can optionally offload the matrix–vector products to GPU hardware.

The main entry point is `solve_sparse_system()` in
`solvers/mod_sparse.f90`.

---

## Direct Solvers

When `gmres = .false.` in the input file, the system is passed directly to one
of the supported sparse-direct libraries.  The library is selected at compile
time via flags (`USE_MUMPS`, `USE_PASTIX`, `USE_STRUMPACK`) and at
runtime by the corresponding `use_*` input flags.


> <b><font color="red">Hint</font></b>:
> Factorization of large sparse matrices is an expensive operation and will create a lot of fill-in resulting in increased memory requirements. Thus, for most cases use of the **iterative solver** is recommended! 

---

## Iterative Solvers

When `gmres = .true.`, a Krylov iterative method is used, preconditioned by the
**block harmonic preconditioner** described in the next section.

### GMRES

Restarted GMRES is implemented in `solvers/mod_gmres2.f90`.

#### Preconditioning side

The preconditioner $M$ can be applied from either side, selected at run time
with `gmres_right_prec`:

| `gmres_right_prec` | System solved | Residual minimised and tested |
|---|---|---|
| `.true.` (default) | $A M^{-1} y = b$, with $x = M^{-1} y$ | $\lVert b - Ax \rVert$, the **true** residual |
| `.false.` | $M^{-1} A x = M^{-1} b$ | $\lVert M^{-1}(b - Ax) \rVert$, the preconditioned residual |

Right preconditioning is the default because GMRES then minimises, and tests
for convergence, the true residual $r = b - Ax$ in the plain Euclidean norm.
Under left preconditioning it works with $M^{-1}r$ instead. Since
$\|r\| \le \|M\|\,\|M^{-1}r\|$, a small preconditioned residual only implies a
small true residual when $\|M\|$ is moderate, and the weighting by $M^{-1}$
changes whenever $M$ does. In JOREK $M$ changes often: the factorisation is
deliberately reused over several time steps (see
[Factorisation Reuse](#factorisation-reuse)), so its quality varies from step
to step.

The choice of side does not make GMRES faster in principle.
$AM^{-1} = M\,(M^{-1}A)\,M^{-1}$, so both preconditioned matrices have the
same eigenvalues; iteration counts differ only through non-normality and the
different norm being minimised.

#### Algorithm

The Arnoldi vectors are orthogonalised by Gram–Schmidt. The variant is chosen
by logical flags hard-coded at the top of `gmres2_driver` (`GSC`, `GSM`,
`GSCI`, `GSMI`). The default, `GSCI`, is classical Gram–Schmidt with
re-orthogonalisation: the projection is repeated, up to three times, as long
as a pass shrinks the vector's norm by more than a factor of two. Givens
rotations keep the upper Hessenberg matrix in triangular form in place.

If the new Arnoldi vector almost vanishes after orthogonalisation, i.e.
$h_{j+1,j} \le 10^{-14}\,\|w_j\|$ where $w_j$ is the vector before
orthogonalisation, the Krylov subspace is numerically invariant and the
restart cycle ends (a *breakdown*). Normally this is a *happy* breakdown: the
small least-squares problem is then solved exactly, the residual estimate
becomes zero, and the solve is reported as converged. It is reported as a
failure only if the entire rotated column vanishes, which makes the
triangular factor singular.

At the end of each cycle of $k$ iterations the solution is updated as

$$x \leftarrow x + M^{-1} V_k\, y_k \quad \text{(right)}, \qquad
x \leftarrow x + V_k\, y_k \quad \text{(left)},$$

where $V_k = [v_1, \dots, v_k]$ holds the Arnoldi vectors and $y_k$ minimises
the residual over the current Krylov subspace. For right preconditioning the
vectors $z_j = M^{-1} v_j$ are stored during the cycle, so the update needs no
further preconditioner solve. This costs `gmres_m` extra stored vectors,
roughly doubling the Krylov memory (the basis itself holds `gmres_m`+1).
Because the update uses the stored $z_j$ instead of applying $M^{-1}$ again,
this is the *flexible* GMRES form: it remains correct even if the
preconditioner changes from one iteration to the next.

#### Convergence

The iteration does not start from zero. Before calling GMRES,
`solve_sparse_system` applies the preconditioner to the right-hand side, so
the initial guess is $x_0 = M^{-1} b$, with residual $r_0 = b - A x_0$.

Let $\rho_k$ be the residual norm that GMRES tracks, and $\beta$ the matching
norm of the right-hand side:

| | $\rho_k$ | $\beta$ |
|---|---|---|
| right | $\lVert b - A x_k \rVert$ | $\lVert b \rVert$ |
| left | $\lVert M^{-1}(b - A x_k) \rVert$ | $\lVert M^{-1} b \rVert$ |

Convergence is declared when

$$\rho_k \le \max\left(\texttt{gmres_tol}\cdot\rho_0,\ \varepsilon\,\beta\right),
\qquad \varepsilon = 10^{-14}.$$

The target is fixed at the first iteration and is not reset by restarts.
Within a cycle, $\rho_k$ is the estimate given by the Givens rotations, which
equals the true value in exact arithmetic; at every restart it is recomputed
from $b - Ax$.

The floor stops a very small `gmres_tol` from iterating on round-off until
`gmres_max_iter`. It is a heuristic, not a guarantee: the smallest residual
GMRES can reach in double precision is of order
$u\,(\|A\|\,\|x\| + \|b\|)$, $u \approx 1.1\times10^{-16}$, which for an
ill-conditioned system can exceed $\varepsilon\,\|b\|$.

`gmres2_driver` returns a `converged` flag, which `solvers/mod_sparse.f90`
passes on as `solver%step_success`. The flag is false when `gmres_max_iter`
iterations are reached without meeting the target, or on a breakdown that is
not a happy one. Each iteration prints $\rho_k$, $\rho_k/\rho_0$
(`res/res0`) and $\rho_k/\beta$ (`res/rhs`).

Key input parameters:

| Parameter | Default | Description |
|---|---|---|
| `gmres_m` | 20 | Restart length: Arnoldi vectors per cycle (capped at `gmres_max_iter`) |
| `gmres_max_iter` | 200 | Maximum **total** number of iterations, counted across restarts |
| `gmres_tol` | 1.d-8 | Tolerance on $\rho_k$ relative to $\rho_0$ (see the floor above) |
| `gmres_right_prec` | .true. | Right (`.true.`) or left (`.false.`) preconditioning |

### BiCGSTAB

BiCGSTAB (in `solvers/mod_bicgstab.f90`) is an
alternative that requires no restart and uses $O(n)$ memory independent of the
iteration count.  It applies the preconditioner **twice per iteration** (once
for the search direction, once for the stabilizer), which is more expensive per
step than GMRES but avoids restarting costs.  Enable it at compile time with the
`USE_BICGSTAB` preprocessor flag.

> <b><font color="red">Hint:</font></b> Use of BiCGSTAB is generally not recommended.

---

## Preconditioner

**Preconditioning** is essential in solving linear systems using iterative methods. The matrices resulting from the MHD formulations in JOREK are highly ill-conditioned and iterative methods do not converge without a proper preconditioner.  
Instead of $Ax = b$, GMRES solves one of two preconditioned systems, the right-preconditioned

$$A M^{-1} y = b, \qquad x = M^{-1} y,$$

or the left-preconditioned

$$M^{-1}A x = M^{-1}b.$$

$M$ is chosen so that systems with $M$ are cheap to solve, and so that $AM^{-1}$ (or $M^{-1}A$) is close to the identity. `gmres_right_prec` selects the side; see [GMRES](#gmres) above.

There exist many "standard" techniques that can be applied to general systems without considering the original problem. These, of course, have the benefit of being "plug-and-play" and work reasonably well for a wide range of problems. When going into very ill-conditioned problems we find most of these approaches to have little to no effect. These cases require more advanced study to create tailored preconditioners to obtain good convergence.

The _block harmonic_ preconditioner in JOREK exploits the **toroidal Fourier structure** of
the MHD system to construct a preconditioner $M$ whose application requires only
independent sparse-direct solves on smaller sub-systems, one per **mode family**.

The preconditioner is set up and applied by
`solvers/mod_preconditioner.f90`.

### Concept

For a problem linearised about an axisymmetric equilibrium, toroidal
harmonics decouple completely: the mode-$n$ rows of the system matrix have no
coupling to mode-$m \ne n$ columns. Nonlinear terms and non-axisymmetric
geometry (e.g. stellarators) do couple harmonics, but in many cases the
dominant coupling is still within a harmonic.  The
preconditioner $M$ is therefore built by **grouping the toroidal modes into
families** and assembling one block-diagonal preconditioner matrix per family
that retains all intra-family coupling while discarding the inter-family terms.
Each mode family is then factorised independently by the selected direct-solver
library.

### Mode Family Distribution

Toroidal modes $0, 1, 2, 3, \ldots$ are partitioned into mode families.
Two strategies are available:

- **Automatic** (`autodistribute_modes = .true.`): mode 0 (axisymmetric) forms
  its own family; remaining modes are paired as $(n_1, n_2), (n_3, n_4), \ldots$
- **Manual**: the user specifies `n_mode_families`, `modes_per_family(:)`, and
  `mode_families_modes(:,:)` explicitly.

#### Overlapping Mode Families

Manual families do not have to form a disjoint partition: the same toroidal
mode can appear in more than one family. Couplings between all modes belonging
to a family are retained in that family's preconditioner block, which can
improve GMRES convergence when important inter-harmonic couplings would
otherwise be discarded.

When families overlap, their solutions contain contributions to the same rows
of the global solution vector. JOREK combines these contributions using the
factor assigned to each family in `weights_per_family`. The weights should be
chosen consistently with the overlap; in the example below every mode belongs
to two families, so every family has weight $0.5$:

```fortran
autodistribute_modes = .false.
n_mode_families = 5

modes_per_family = 1, 3, 4, 4, 2
mode_families_modes(1,:) = 1
mode_families_modes(2,:) = 1, 2, 3
mode_families_modes(3,:) = 2, 3, 4, 5
mode_families_modes(4,:) = 4, 5, 6, 7
mode_families_modes(5,:) = 6, 7

weights_per_family = 0.5, 0.5, 0.5, 0.5, 0.5
```

The values in `mode_families_modes` are the one-based `i_tor` mode indices used
internally by JOREK. Increasing the number of modes per family retains more of
the original matrix coupling, but also produces larger and more expensive
blocks to factorise.

MPI ranks are distributed among families in the same way:

- **Automatic** (`autodistribute_ranks = .true.`): ranks are distributed as
  equally as possible across families.
- **Manual**: specified via `ranks_per_family(:)`.

Because families can have very different block sizes, an equal rank
distribution may give poor load balance. A useful starting heuristic is to
assign ranks approximately in proportion to the square of the number of modes
in each family. For the five overlapping families above, one possible setup is:

```fortran
autodistribute_ranks = .false.
ranks_per_family = 2, 8, 16, 16, 4
```

The total number of MPI ranks used for the run must equal the sum of
`ranks_per_family`. The optimal distribution depends on the matrix structure,
the direct solver, and the machine, so this initial allocation may require
tuning.

For further details, see [Holod et al., *Enhanced preconditioner for JOREK MHD
solver*, Plasma Physics and Controlled Fusion **63**, 114002
(2021)](https://iopscience.iop.org/article/10.1088/1361-6587/ac206b).

<div style="display: flex; justify-content: space-evenly; align-items: flex-start; text-align: center;">
  
  <figure style="width: 30%; margin: 0;">
    <img src="assets/solver/pc_matrix.png" alt="Single block harmonic" style="width: 100%;">
    <figcaption><i>Single harmonic pair per block</i></figcaption>
  </figure>

  <figure style="width: 31.5%; margin: 0;">
    <img src="assets/solver/pc_matrix_coupled.png" alt="Grouped block harmonic" style="width: 100%;">
    <figcaption><i>Grouped harmonic families</i></figcaption>
  </figure>

</div>

### Preconditioner Solve Workflow

Each GMRES / BiCGSTAB preconditioner application performs the following steps:

1. **Scatter RHS** — extract the rows belonging to this mode family from the
   global residual vector (via pre-computed `row_index` mapping).
2. **Solve** — the family's direct solver (MUMPS / PaStiX / STRUMPACK)
   factorises or reuses, and solves the preconditioner matrix for the local
   RHS.
3. **Gather solution** — contributions from all families are reduced by
   `MPI_AllReduce` (sum) into the global solution vector; each family's rows are
   weighted by `row_factor` (normally 1, or set from `weights_per_family` for
   manually defined families).

### Factorisation Reuse

Refactoring the preconditioner at every time step is expensive. Thus, we can reuse the factorized preconditioner matrices over multiple timesteps, where they still provide a good approximation of the linear system.  After the first time step, JOREK reuses the existing factorisation (the
`solve_only` path) when both

$$\texttt{iter_gmres} + \texttt{iter_prev} \le 2 \times \texttt{iter_precon}
\quad \text{and} \quad
\texttt{n_since_update} < \texttt{max_steps_noUpdate}$$

hold, i.e. the last two solves averaged at most `iter_precon` GMRES iterations
and the factorisation is younger than `max_steps_noUpdate` steps (both are
input parameters). If either condition fails, the preconditioner matrix is
reassembled and refactorised. Within a Newton loop the preconditioner is never
updated after the first Newton iteration.

### Preconditioner Matrix Assembly

Two assembly strategies are available (selected by compile-time flag
`DIRECT_CONSTRUCTION`):

- **Direct construction**: each MPI rank assembles its part of the
  preconditioner matrix directly at the finite-element level, without the
  distribution step.
- **Distributed construction**: the full system matrix is communicated via
  `MPI_ALLTOALLV` to the rank groups that own each mode family, which then
  extract their rows locally.

---

## PETSc Integration

> **Note:** PETSc-based solver paths are under active development.
> This section will be expanded once the implementation stabilises.

When compiled with `USE_PETSC`, an alternative solver path is
available that uses PETSc's KSP framework for both direct and iterative solves.
The PETSc integration layer lives in `solvers/mod_petsc.f90` and exposes:

- A persistent `KSP` context that is reused across time steps.
- An MPIBAIJ system matrix that can be filled from JOREK's native BCSR data.
- Plug-in preconditioners implemented as PETSc `PCSHELL` or `PCFIELDSPLIT`
  objects.

Details of the PETSc-based preconditioners and their configuration will be
documented separately.

---

## Further reading

#### JOREK solver and numerical methods
- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM (2003). [Freely available from the author.](https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf)
- T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM (2006).
- I. Holod *et al.*, "Enhanced preconditioner for JOREK MHD solver," arXiv:2101.08646 (2021).
- A. Quinlan, V. Dwarka, I. Holod, M. Hoelzl, "Towards Robust Solvers for Nuclear Fusion Simulations Using JOREK: A Numerical Analysis Perspective," arXiv:2308.16124 (2023).

#### Direct Solver back-ends
- **MUMPS** — [mumps-solver.org](https://mumps-solver.org/)
- **PaStiX** — [solverstack.gitlabpages.inria.fr/pastix](https://solverstack.gitlabpages.inria.fr/pastix/)
- **STRUMPACK** — [portal.nersc.gov/project/sparse/strumpack](https://portal.nersc.gov/project/sparse/strumpack/)
- **PETSc** — [petsc.org](https://petsc.org/)
