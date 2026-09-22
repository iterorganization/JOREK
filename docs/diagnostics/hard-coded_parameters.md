---
title: "Hard-Coded Parameters in JOREK"
nav_order: 16
render_with_liquid: false
parent: "Installation & Compilation"
---

## Physics Model

The first "hard-coded parameter" is the **JOREK physics model number**. It is set in the Makefile.inc file.

The other hard-coded parameters are set by the mod_parameters.f90 file located in the model directory models/modelXXX/:

| Parameter | Description |
| --- | --- |
| n_tor | Number of toroidal mode numbers (cos+sin counted separately) |
| n_period | Periodicity (n_period=1 means full torus) |
| n_plane | Number of toroidal planes; n_plane must be >= 2*(n_tor-1) |
| n_vertex_max | do not change |
| n_nodes_max | Maximum number of grid nodes |
| n_elements_max | Maximum number of grid elements (2D Bezier) |
| n_boundary_max | Maximum number of boundary elements (1D Bezier) |
| n_pieces_max | Maximum number of pieces in a flux surface representation |

When changing hard-coded parameters, you should always execute:
  make clean

When changing model you should execute:
  make cleanall

## Script config.sh

The script
  ./util/config.sh
called from the trunk allows to **show the values of the hard-coded parameters**. You will get an output like this:
  =======================
  Makefile.inc:
  -----------------------
    model303
    n_tor = 1
    n_period = 1
    n_plane = 8
    n_vertex_max = 4
    n_nodes_max = 60001
    n_elements_max = 60001
    n_boundary_max = 1001
    n_pieces_max = 6001
  =======================

The same script can also **modify the hard-coded parameters**:
  ./util/config.sh model=303 n_tor=3 n_period=8 n_plane=8

**Important:** When changing the model with the script, make sure that the other parameters are correct as well, since they are set separately for each model in the model-specific mod_parameters.f90 files.
