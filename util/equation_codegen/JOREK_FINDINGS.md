# Model-600 Jacobian: missing and inconsistent terms

This file records the differences between the volume terms assembled by
`models/model600/mod_elt_matrix_fft.f90` and the automatically linearized
weak forms in `src/jorek_equations/model600.py`.

Every finding below was obtained by generating the two term reports and
comparing them block by block:

```bash
python3 examples/export_model600_terms.py
python3 examples/diff_model600_reports.py
```

`diff_model600_reports.py` compares each `rhs_ij`/`amat` assignment as a
multiset of monomials and prints the algebraic residual, so a block that only
reorders its terms is reported as identical.  The residual is evaluated after
rewriting every `_s`/`_t` derivative through the chain rule, because the
element routine spells the same poloidal bracket sometimes as
`a_s*b_t - a_t*b_s` and sometimes as `xjac*(a_x*b_y - a_y*b_x)` — even between
a residual and its own tangent.  A block whose two sides differ only by that
rewrite is reported as `SAME`.

The original audit reported seventeen findings. Findings 3, 5–8 and 10–17
have since been fixed in the Fortran (one commit per finding, `fix
linearization: finding N, ...`), and their sections have been removed from
this file. The neutral density row `rhon` was added afterwards and brought
findings 18 and 19, which have also since been fixed. What is left disagrees
in 24 blocks, for two reasons:

- **Open: the inward pinch (findings 4 and 9).** These are real
  inconsistencies in the Jacobian, but the pinch implementation as a whole
  needs revisiting, so they are left to its developer rather than patched
  here. They are described below.
- **Accepted: tauIC with impurities (findings 1 and 2).** These terms only
  exist when the diamagnetic terms (`tauIC /= 0`) and impurities are both
  active, and that combination is not supported. They are kept in the
  reference as known differences, not as bugs to fix. See the section below.

Status of the exported rows:

| Row | Residual | Jacobian |
|---|---|---|
| `psi` | reproduced | finding 2 (tauIC × impurities only) |
| `u` | reproduced | finding 1 (tauIC × impurities only) |
| `zj` | reproduced | reproduced |
| `w` | reproduced | reproduced |
| `rho` | reproduced | findings 1 (tauIC × impurities only) and 4 |
| `vpar` | reproduced | findings 4 and 9 |
| `rhoimp` | reproduced | reproduced |
| `rhon` | reproduced | reproduced |
| `Ti` | reproduced | reproduced |
| `Te` | reproduced | reproduced |
| `T` | reproduced | reproduced |

Every residual is now reproduced exactly, so none of the remaining differences
changes the converged solution. They only affect the Newton tangent.

The NEO branch is excluded from the default reports (`--include-neo` enables
it) and has not been audited here.

---


## Finding 4 — the density inward pinch is not differentiated with respect to psi

**Where:** the `var_rho` block, `amat(var_rho,var_psi)`.

**What happens:** the density residual contains the weak form of
`-div(r0 * V_pinch)` with `V_pinch = -V_prof_pinch * grad(psi)/|grad(psi)|`:

```fortran
- V_prof_pinch / sqrt(psi_grad2) * (v_x * ps0_x + v_y * ps0_y) * r0 * BigR * xjac * tstep * factor(var_rho,14)
```

with `psi_grad2 = ps0_x**2 + ps0_y**2`. The pinch direction is therefore a
function of the evolved poloidal flux, but `amat(var_rho,var_psi)` contains no
corresponding tangent. Every other flux-dependent quantity of the same
equation *is* differentiated there — `BB2` through `BB2_psi`, the parallel
gradients through `Bgrad_rho_star_psi`, `Bgrad_rho_psi` and
`Bgrad_rhoimp_psi` — which is what makes the omission stand out.

**Suggested fix:** the missing contribution is the transverse projection of
the trial flux gradient,

```fortran
+ V_prof_pinch / sqrt(psi_grad2) * r0 * BigR * xjac * theta * tstep                       &
    * ( (v_x * psi_x + v_y * psi_y)                                                       &
      - (ps0_x*v_x + ps0_y*v_y) * (ps0_x*psi_x + ps0_y*psi_y) / psi_grad2 )
```

(the first bracket from varying `grad(psi)` in the dot product, the second from
varying `1/|grad(psi)|`).

The parallel-velocity equation advects `v_par` with the same pinch velocity
and has the same gap: `amat(var_vpar,var_psi)` carries no pinch tangent
either.

**Scale:** 6 monomials in `amat(var_rho,var_psi)` and 6 in
`amat(var_vpar,var_psi)`, in each temperature branch.

**Note:** `V_prof_pinch`, `D_prof`, `D_prof_imp`, `D_par_local`,
`D_par_local_imp`, `tau_sc` and `D_perp_num_psin` are profile work values that
also depend on the state through `psi_norm` (and, for `tau_sc`, through the
pressure). The element routine freezes all of them, and the generator follows
that convention, so they produce no report differences. The pinch term is
different: its flux dependence is explicit in the assembled expression rather
than hidden inside a work value.

---

## Finding 9 — the parallel-velocity inward-pinch tangent has the wrong sign

**Where:** `amat(var_vpar,var_rho)` and `amat(var_vpar,var_vpar)`.

**What happens:** the residual contribution is

```fortran
rhs_ij(var_vpar) = rhs_ij(var_vpar) &
    + V_prof_pinch / sqrt(psi_grad2) * (ps0_x*vpar0_x + ps0_y*vpar0_y) &
            * r0 * v * BigR * xjac * tstep * factor(var_vpar,13)
```

and its two tangents repeat that sign:

```fortran
amat(var_vpar,var_rho)  = ... + V_prof_pinch / sqrt(psi_grad2) * (ps0_x*vpar0_x + ps0_y*vpar0_y) * rho  * v * BigR * xjac * theta * tstep
amat(var_vpar,var_vpar) = ... + V_prof_pinch / sqrt(psi_grad2) * (ps0_x*vpar_x  + ps0_y*vpar_y ) * r0   * v * BigR * xjac * theta * tstep
```

The JOREK convention is `AMAT = -theta*tstep*dB`, so both should be negative.
The density equation gets this right for the same pinch term: its residual
contribution is negative and `amat(var_rho,var_rho)` is positive.

**Suggested fix:** flip the sign of both `amat(var_vpar,...)` pinch terms.

**Note:** the residual sign itself also deserves a look. With
`V_pinch = -V_prof_pinch * grad(psi)/|grad(psi)|` as documented in the source
comment, the advection `rho * V_pinch . grad(v_par)` is
`-V_prof_pinch/|grad(psi)| * rho * grad(psi).grad(v_par)`, i.e. negative,
while the assembled term is positive. The generator follows the assembled
residual, so this note is a reading of the source comment rather than a
report difference.

**Scale:** 2 monomials in each of `amat(var_vpar,var_rho)` and
`amat(var_vpar,var_vpar)`, in each temperature branch.

---

## Accepted differences: tauIC with impurities (findings 1 and 2)

The diamagnetic terms scale with `tauIC`. Their tangents are written for the
main-ion and electron pressures `r0*Ti0` and `r0*Te0`, while the residuals use
the full pressures `(r0 + rimp0*alpha_i)*Ti0` and `(r0 + rimp0*alpha_e)*Te0`.
The generator therefore produces impurity-pressure terms that the element
routine does not have:

- **Finding 1**: the impurity ion pressure in the diamagnetic terms of the
  momentum and density equations: `amat(var_u,var_Ti|T|rhoimp)` and
  `amat(var_rho,var_Ti|T|rhoimp)`.
- **Finding 2**: the impurity electron pressure in the diamagnetic coupling
  of the induction equation: `amat(var_psi,var_Te|T|rhoimp)`,
  `amat_n(var_psi,var_Te|T)`, and an `amat_n(var_psi,var_rhoimp)` that the
  element routine does not assemble at all.

Every one of these monomials carries both `tauIC` and an impurity factor
(`rhoimp`, `rimp0` or its derivatives). Model 600 does not support running
with impurities and `tauIC /= 0` at the same time, so the missing terms never
contribute in a supported configuration. They are recorded in the reference as
accepted differences, not bugs. The combination should be rejected at input
time, but nothing in model 600 does this yet. If it ever becomes supported,
these tangents must be completed first; the full derivation is in the history
of this file (commit `f4f8a4e98`).

---

## Complete inventory

Every report line that is blank on one side belongs to one of the four
findings above. There is nothing unexplained left.

| Finding | Report lines | Blocks (each in both temperature branches) |
|---:|---:|---|
| 1 — impurity ion pressure, tauIC terms (accepted) | 96 | `amat(var_u,var_rhoimp)` 44, `amat(var_u,var_t)` 22, `amat(var_u,var_ti)` 22, `amat(var_rho,var_rhoimp)` 4, `amat(var_rho,var_t)` 2, `amat(var_rho,var_ti)` 2 |
| 2 — impurity electron pressure, tauIC terms (accepted) | 30 | `amat(var_psi,var_rhoimp)` 10, `amat(var_psi,var_t)` 8, `amat(var_psi,var_te)` 8, `amat_n(var_psi,var_rhoimp)` 2, `amat_n(var_psi,var_t)` 1, `amat_n(var_psi,var_te)` 1 |
| 4 — pinch not differentiated with respect to psi (open) | 24 | `amat(var_rho,var_psi)` 12, `amat(var_vpar,var_psi)` 12 |
| 9 — sign of the vpar pinch tangent (open) | 16 | `amat(var_vpar,var_rho)` 8, `amat(var_vpar,var_vpar)` 8 |
| **total** | **166** | 24 blocks: 8 source-only and 158 generated-only lines |

The 8 source-only lines are the wrong-sign pinch terms of finding 9, two per
block in both temperature branches. Each has a generated partner with the
opposite sign. All other lines are terms the generator produces and the
element routine does not.

---

## Regression benchmark

The inventory above is frozen in
[`reference/model600_discrepancies.json`](reference/model600_discrepancies.json)
and checked by [`final_test.py`](final_test.py) (or `run_test.sh`, which also
sets up `sympy`). Each block of the reference names the finding it belongs
to, so a failure points straight back at this document:

```bash
./run_test.sh              # or: ./final_test.py
```

Fixing finding 4 or 9 in the Fortran will make the benchmark fail with `-`
lines for the monomials that stopped disagreeing. That is the intended signal:
re-run with `--update`, then remove the finding from this file.

## Reproducing

```bash
cd util/equation_codegen
module load sympy/1.14.0-gfbf-2025b        # or any Python with sympy
python3 examples/export_model600_terms.py
python3 examples/diff_model600_reports.py
```

The **residual** printed under a `DIFF` block is the authoritative difference.
It is computed in the element basis, where the two spellings of a poloidal
bracket coincide. It is displayed with `f_s -> f_x`, `f_t -> f_y`,
`xjac -> 1` so that it reads as a physical expression.

The tool also verifies its own alignment. The blank cells of every block must
account for that block's multiset difference, and any surplus is reported as
`MISALIGNED`. In the aligned Markdown reports, a blank cell on the generated
side is therefore a genuinely missing term.
