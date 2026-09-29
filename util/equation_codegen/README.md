# JOREK equation code generator

This package checks the hand-written linearization in a JOREK element routine
against one derived symbolically from the same weak form. It does four things:

1. **Read** the assembled `rhs_ij`/`amat` terms out of a model's
   `mod_elt_matrix_fft.f90`.
2. **Generate** the residual and the Jacobian blocks from the weak form, by
   directional differentiation of the equations written in the DSL.
3. **Compare** the two, monomial by monomial.
4. **Report** the comparison, so that a mismatch is readable rather than a
   wall of algebra.

The symbolic layer is described in [`SPECIFICATION.md`](SPECIFICATION.md).
Equations are written with bare fields (`psi`, `u`, `rho`, ...); the
linearization stage decides when those become background values and when the
differentiated one becomes a trial function.

## Layout

| Path | Role |
|---|---|
| `src/jorek_equations/symbols.py`, `operators.py`, `external.py` | fields, roles, spatial operators, externally supplied quantities |
| `src/jorek_equations/equations.py`, `linearization.py`, `channels.py` | weak equations, the Gateaux derivative, the FFT channel split |
| `src/jorek_equations/fortran.py` | printing a symbolic term with JOREK's names |
| `src/jorek_equations/fortran_source.py` | step 1 and step 3: read the element routine, normalize its work variables, compare |
| `src/jorek_equations/model199.py`, `model600.py` | the equations themselves |
| `src/jorek_equations/model600_markdown.py` | step 4: the aligned Markdown reports |
| `final_test.py`, `reference/` | the regression benchmark |

## Setup and tests

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -m unittest discover -s tests -v
```

## Model 199

```bash
.venv/bin/python examples/check_model199.py
```

The comparison is symbolic rather than textual. A mismatch raises an error
that includes the source expression, the generated expression, and their
algebraic difference.

## Model 600 reports

`export_model600_terms.py` writes the two line-aligned Markdown reports for
the `psi`, `u`, `zj`, `w`, `rho`, `vpar`, `rhoimp`, `Ti`, `Te` and `T` rows:

```bash
.venv/bin/python examples/export_model600_terms.py
meld reports/model600_fortran_terms.md reports/model600_generated_terms.md
```

Repeat `--equation` to export fewer rows while developing one of them:

```bash
.venv/bin/python examples/export_model600_terms.py --equation psi
```

`Ti` and `Te` exist only in the two-temperature branch of the element routine
and `T` only in the other one, so each is exported into its own section.

To compare the two reports block by block instead of line by line:

```bash
.venv/bin/python examples/diff_model600_reports.py
.venv/bin/python examples/diff_model600_reports.py --quiet   # verdicts only
```

Each assignment block is compared as a multiset of monomials, and a block that
differs is printed together with the algebraic residual `source-generated`.
Reordering alone therefore does not show up as a difference. The residual is
taken after rewriting every `_s`/`_t` derivative through the chain rule, so a
poloidal bracket written `a_s*b_t - a_t*b_s` in one report and
`xjac*(a_x*b_y - a_y*b_x)` in the other cancels; such a block is reported as
`SAME` rather than `DIFF`.

## Temperature conventions

The element routine builds every pressure in `construct_pressure`, which is
written once for both temperature models and always uses the species
temperatures `Ti0`/`Te0` and the per-species impurity coefficients `alpha_i`,
`alpha_e` and `alpha_e_bis`. The one-temperature model evolves the *total*
temperature and enters that routine with

```
Ti0 = Te0 = T0/2
```

so a single-species temperature is half of the evolved field, and the stored
one-temperature closure values are the means

```
alpha_imp     = (alpha_i + alpha_e)/2
alpha_imp_bis = (alpha_i + alpha_e_bis)/2
alpha_imp_tri = alpha_e_tri/4
```

The comparison therefore expands both reports in the two-species basis. The
DSL follows the same convention: `_diamagnetic_pressure` and the induction
equation's electron pressure use `T/2` in the one-temperature branch. This
replaced two earlier "legacy tangent" overrides that reproduced JOREK's
factor of one half by hand; that half is physical, not a legacy quirk.

## The `final_test` benchmark

The discrepancies between the element routine and the generated linearization
are frozen in
[`reference/model600_discrepancies.json`](reference/model600_discrepancies.json).
For every assignment that disagrees, it records the monomials the element
routine has and the generator does not, the monomials the generator has and
the element routine does not, and which findings of
[`JOREK_FINDINGS.md`](JOREK_FINDINGS.md) that block belongs to. Those two
lists are the complete statement of the disagreement: everything else in the
reports is identical on both sides.

```bash
./final_test.py                    # recompute and compare; exit 0 on a match
./final_test.py --equation vpar    # only these rows, for a quick check
./final_test.py --update           # rewrite the reference after an intended change
```

The benchmark exports the reports into a temporary directory, so it never
touches `reports/`. It takes a few minutes, which is why it is a standalone
script rather than part of `python -m unittest discover`; the unit tests cover
its comparison logic and the integrity of the reference file.

A failure prints one line per changed monomial: `-` for a discrepancy the
reference expects that the run no longer produces, `+` for a new one. Both
directions matter. A `-` line means either that somebody fixed the Fortran, in
which case update the reference and strike the finding from
`JOREK_FINDINGS.md`, or that the generated equation drifted and stopped
reproducing a term it used to reproduce. A `+` line means a new disagreement:
a change in the element routine, or a regression in the DSL. The script also
warns when `mod_elt_matrix_fft.f90` no longer has the SHA-256 recorded in the
reference, which usually explains the rest of the output.

## Known JOREK differences

With the conventions above, all ten rows — `psi`, `u`, `zj`, `w`, `rho`,
`vpar`, `rhoimp`, `Ti`, `Te` and `T` — are reproduced by the generator, and
174 of the 233 exported assignment blocks are identical. The remaining 59 are
terms the generator produces and the element
routine does not, or coefficients that disagree. They are recorded, with their
suggested Fortran fixes, in [`JOREK_FINDINGS.md`](JOREK_FINDINGS.md):

1. the impurity part of the ion pressure is not differentiated in the
   diamagnetic terms of the momentum and density equations
   (`amat(var_u,var_Ti)`, `amat(var_u,var_T)`, `amat(var_u,var_rhoimp)` and
   the same three columns of `var_rho`);
2. the impurity part of the electron pressure is not differentiated in the
   induction equation (`amat(var_psi,var_Te)`, `amat_n(var_psi,var_Te)`,
   `amat(var_psi,var_rhoimp)`, and an `amat_n(var_psi,var_rhoimp)` that the
   element routine does not have at all);
3. the one-temperature `amat(var_u,var_T)` diamagnetic-viscosity tangent is a
   factor two too large;
4. the density inward-pinch term is not differentiated with respect to `psi`
   (`amat(var_rho,var_psi)`);
5. the temperature dependence of `alpha_e` is not differentiated in the
   density ionization/recombination sources (`amat(var_rho,var_Te)`,
   `amat(var_rho,var_T)`), although the momentum equation does differentiate
   it in the same sources;
6. the parallel-velocity time term is linearized inconsistently: the `vpar`
   column and the history keep only the toroidal part of `B**2`, and the
   `psi` column carries half of its variation;
7. `BB2` is not differentiated in the `tgnum_vpar` tangent
   (`amat(var_vpar,var_psi)`, `amat_k(var_vpar,var_psi)`);
8. the toroidal channel of the parallel-parallel viscosity tangent is missing
   (`amat_n(var_vpar,var_vpar)`, `amat_kn(var_vpar,var_vpar)`);
9. the parallel-velocity inward-pinch tangent repeats the sign of its own
   residual instead of flipping it;
10. the impurity parallel diffusivity loses its shock-capturing part
    `D_par_imp_sc_num*tau_sc` in the toroidal channel of the residual
    (`rhs_ij_k(var_rhoimp)`) — the only finding that affects the converged
    solution rather than only the Newton tangent;
11. the `tgnum_Ti` poloidal-velocity tangent carries `BigR**2` where its own
    residual carries `BigR**3`, in all four of its columns;
12. three smaller omissions in the ion energy tangents: the `BB2` flux
    dependence of the kinetic-coupling term, the density derivative of the
    perpendicular conductivity in the toroidal channel, and the electron
    temperature dependence of the recombination rate;
13. four ionization-energy tangent lines of `amat(var_Te,var_Te)` and
    `amat_k(var_Te,var_Te)` are missing their `theta` factor, so those entries
    are `1/theta` times too large;
14. `amat(var_Te,var_Te)` uses `alpha_e` where `amat_k(var_Te,var_Te)` uses
    `alpha_e_bis` for the same quantity, two lines apart;
15. four omissions in the electron and total energy tangents: the temperature
    dependence of six ionization-energy terms, uncorrected densities in
    `amat(var_Te,var_rhon)`, three radiation/ionization derivatives in
    `amat(var_Te,var_rhoimp)`, and `alpha_e` in the friction sources;
16. `rhs_ij(var_T)` convects the pressure with the *neutral* density `rn0`
    where every analogous line uses `rimp0` — the second finding that affects
    the converged solution;
17. `amat_n(var_T,var_vpar)` drops the impurity part of the pressure that both
    two-temperature equations keep.

A source monomial and its generated counterpart are aligned on the same report
line even when only their numeric coefficient or their power of `BigR` differs, so a term the element
routine scales differently shows up side by side in Meld rather than as a
source line with an empty generated cell followed by an unmatched generated
line at the end of the block.

The element routine is not consistent about how it spells a poloidal bracket:
`rhs_ij(var_vpar)` writes the parallel kinetic-energy flux as
`a_s*b_t - a_t*b_s` while `amat(var_vpar,var_psi)` writes the very same group
as `xjac*(a_x*b_y - a_y*b_x)`, and the two expand into different monomials.
The parallel-velocity equation is therefore generated in both spellings and
each assignment is aligned against whichever one its own source block uses.
A blank cell on the generated side of the report is then a genuinely missing
term rather than a change of coordinates.

Matching runs in three passes over a whole block — exact monomials first, then
the ones that differ only by a coefficient, then the ones that also differ by
a power of `BigR`. A single interleaved pass lets an early source monomial
with no exact partner consume, through one of the relaxed keys, a generated
monomial that a later source monomial matches exactly; the pair then drifts
apart even though the block agrees. A fourth pass pairs terms that differ by the implicitness factor `theta` as
well, since that is a scalar scaling like any other. `diff_model600_reports.py`
checks the result: the blank cells of every block must account for its
multiset difference, and any surplus is reported as `MISALIGNED`.

Assignments and terms follow their order in
`models/model600/mod_elt_matrix_fft.f90`. A source term that is not available
from the equation generator has an empty line at the same position in the
generated report. Generated-only terms are appended to their corresponding
assignment block with an empty source line. NEO AMAT terms remain explicitly
source-only until their symbolic linearization is enabled.
