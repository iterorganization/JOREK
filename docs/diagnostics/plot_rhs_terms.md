---
title: "Plot Equation Terms in VTK"
nav_order: 4
parent: "Diagnostics and post-processing"
layout: default
render_with_liquid: false
---

# Plotting Separate Equation Terms in VTK

This diagnostic is available only for model 600.

## Overview

This `jorek2_postproc` tool reads a restart file and produces a VTK file in which the terms of each equation in `elm_matrix` are separated. It can provide insight into the dominant physical terms and help diagnose numerical problems.

Features include:

- **Code duplication is avoided:** The terms are taken directly from the right-hand side of `mod_elm_matrix_fft`.
- The user can select a single toroidal harmonic with `only_itor`. If it is not specified, all harmonics are summed.
- The user can select a poloidal plane with `vtk_phi_value` and the VTK resolution with `nsub_vtk`.
- The user can calculate the terms of only one equation. If no equation is selected, terms from all equations are produced.
- Terms are labelled according to their physical meaning. The beginning of each label identifies the equation or variable, while the end identifies the particular term in that equation.

**Important:** Terms that are integrated by parts require additional work. Boundary integrals are still missing when `elm_matrix` is called by this diagnostic.

## Usage

This tool is implemented in the [interactive diagnostic tool `jorek2_postproc`](/JOREK/howto/introduction_to_jorek_diagnostics.html#jorek2_postproc). See that page for compilation and general usage instructions.

Unlike most `jorek2_postproc` diagnostics, `RHS_terms_vtk` cannot be run interactively. Create a post-processing script with the following structure:

```text
namelist INPUT_FILE
set nsub_vtk NUMBER_OF_SUBDIVISIONS
set vtk_phi_value TOROIDAL_ANGLE
set only_itor TOROIDAL_HARMONIC
for step FIRST_STEP to LAST_STEP do
  RHS_terms_vtk EQUATION_NUMBER
done
```

The parameters and command arguments are:

- `nsub_vtk`: Number of subdivisions of the JOREK elements used for VTK visualization.
- `vtk_phi_value`: Toroidal angle, in radians, of the poloidal plane to calculate.
- `only_itor`: Selects a single toroidal harmonic. Remove this line to sum all harmonics.
- `EQUATION_NUMBER`: Selects a single equation. Omit this argument from `RHS_terms_vtk` to calculate the terms of every equation.

Run the script using one MPI process:

```bash
export OMP_NUM_THREADS=NUMBER_OF_THREADS
mpirun -n 1 jorek2_postproc < POSTPROCESSING_SCRIPT
```

Do not use more than one MPI process. Large cases may need to be submitted as a batch job, still with only one MPI process. The resulting VTK files are written to the `postproc/` directory.

On one specific Marconi configuration, a problem was resolved by calling `jorek2_postproc` directly without `mpirun`.
