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
| `src/jorek_equations/latex_render.py`, `model600_docs.py`, `model600_notation.py` | the model-600 documentation, generated from `model600.py` |
| `final_test.py`, `reference/` | the regression benchmark |

## Setup and tests

The only dependency is `sympy`. No package install is needed, because the
scripts add `src/` to `sys.path` themselves. On the ITER cluster:

```bash
module load sympy/1.14.0-gfbf-2025b
python3 -m unittest discover -s tests -v
```

Any other Python with `sympy` works too, for example a virtual environment
with `pip install sympy`. `run_test.sh` tries these in order: the system
Python, then `module load`, then a local `.venv`.

## Model 199

```bash
python3 examples/check_model199.py
```

The comparison is symbolic rather than textual. A mismatch raises an error
that includes the source expression, the generated expression, and their
algebraic difference.

## Model 600 documentation

[`docs/physics/base_fluid_models/RMHD/weak_form.md`](../../docs/physics/base_fluid_models/RMHD/weak_form.md)
documents the model-600 equations: every equation written out in full, the
shared operators, and each equation's term groups. Its math is generated from
`model600.py`, which stays the only place the equations are written; the
generator reads the Python source (it does not evaluate it), so each
`_helper` appears as a named operator and each commented block of a `B = (...)`
sum becomes a term group. How each name is written in LaTeX, and the titles and
descriptions of the operators, are in `model600_notation.py`. Only the regions
between `BEGIN GENERATED` and `END GENERATED` are rewritten; the prose around
them is hand-written.

```bash
python3 examples/model600_docs.py render   # rewrite the generated regions
python3 examples/model600_docs.py check    # exit 1 if the page is stale
```

The unit tests and `run_test.sh` run `check`, so a change to `model600.py`
that is not rendered into the page fails the build. A new equation or helper
must also be added to `model600_notation.py`; `render` says so if it is missing.

## Model 600 reports

`export_model600_terms.py` writes the two line-aligned Markdown reports for
the `psi`, `u`, `zj`, `w`, `rho`, `vpar`, `rhoimp`, `rhon`, `Ti`, `Te` and `T` rows:

```bash
python3 examples/export_model600_terms.py
meld reports/model600_fortran_terms.md reports/model600_generated_terms.md
```

Repeat `--equation` to export fewer rows while developing one of them:

```bash
python3 examples/export_model600_terms.py --equation psi
```

`Ti` and `Te` exist only in the two-temperature branch of the element routine
and `T` only in the other one, so each is exported into its own section.

To compare the two reports block by block instead of line by line:

```bash
python3 examples/diff_model600_reports.py
python3 examples/diff_model600_reports.py --quiet   # verdicts only
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

With the conventions above, all eleven rows (`psi`, `u`, `zj`, `w`, `rho`,
`vpar`, `rhoimp`, `rhon`, `Ti`, `Te` and `T`) are reproduced by the
generator. The first audit reported seventeen findings for the first ten
rows. Most have since been fixed in the Fortran, and every residual
(`rhs_ij`) now agrees exactly. The assignment blocks that still differ are all
Jacobian tangents, documented in [`JOREK_FINDINGS.md`](JOREK_FINDINGS.md):

- **Open, left to the pinch developer**:
  - Finding 4: the density and parallel-velocity inward-pinch terms are not
    differentiated with respect to `psi` (`amat(var_rho,var_psi)`,
    `amat(var_vpar,var_psi)`).
  - Finding 9: the parallel-velocity pinch tangent repeats the sign of its
    own residual instead of flipping it (`amat(var_vpar,var_rho)`,
    `amat(var_vpar,var_vpar)`).
- **Accepted, not bugs**:
  - Findings 1 and 2: the impurity parts of the ion and electron pressures
    are not differentiated in the `tauIC` diamagnetic tangents of the `u`,
    `rho` and `psi` equations.
  - Model 600 does not support impurities together with `tauIC /= 0`, so
    these terms never contribute in a supported run.

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
