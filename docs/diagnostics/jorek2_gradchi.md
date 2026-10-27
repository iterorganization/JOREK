---
title: "jorek2_gradchi_domm / jorek2_gradchi_full"
nav_order: 37
render_with_liquid: false
parent: "General"
---

# jorek2_gradchi_domm and jorek2_gradchi_full

These two programs evaluate the vacuum magnetic field of the stellarator models, $\nabla\chi$, at arbitrary $(x,y,z)$ points. They differ in what part of $\chi$ they include.

| Program | Evaluates | Needs |
| :--- | :--- | :--- |
| `jorek2_gradchi_domm` | Dommaschk potentials only | namelist and Dommaschk file, no restart. Build with `USE_DOMM=1` |
| `jorek2_gradchi_full` | full $\chi$ as used in the run (Dommaschk part plus finite-element correction or external field) | namelist and restart file |

`jorek2_gradchi_domm` is fast and needs no restart, so other codes can call it directly. CARIDDI uses it to obtain $\nabla\chi$ on its control surface for the stellarator free-boundary coupling (see [Running CARIDDI](/JOREK/physics/model_extensions/freebound/running_CARIDDI.html)). For tokamak models it returns $F_0\nabla\phi$.

Compile them with the same hard-coded parameters as the JOREK run:

    make -j 8 jorek2_gradchi_domm jorek2_gradchi_full

## Usage

    ./jorek2_gradchi_domm [--xyz xyz.dat] [--out fields_xyz.dat] < namelist
    ./jorek2_gradchi_full [--xyz xyz.dat] [--out fields_xyz.dat] [--restart jorek_restart] < namelist

The values in brackets are the defaults. `--restart` takes the restart name without the `.h5` extension.

## Input

The points in metres, one per line. The first line may hold the number of points. Empty lines and lines starting with `#` are skipped.

    number_of_points
    x_1       y_1       z_1
     :         :         :
    x_n       y_n       z_n

## Output

One header line, then one line per input point in input order

    x  y  z  gradchi_x  gradchi_y  gradchi_z

with the gradient in Tesla. `jorek2_gradchi_full` adds three columns with the element index and the local coordinates `s`, `t` of each point. Points outside the mesh get element index `-1`. Their value is the Dommaschk part alone in a `USE_DOMM` build, where this is exact, and `NaN` otherwise.

Coordinates follow the JOREK convention $x = R\cos\phi$, $y = -R\sin\phi$, $z = Z$.

Errors in the input points, the command line or the files stop the programs with a non-zero exit code. The output file is emptied at start, so a run that stops early never leaves the results of a previous run behind.
