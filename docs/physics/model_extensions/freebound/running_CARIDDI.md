---
title: "Running and setting up CARIDDI"
nav_order: 5
parent: "Free Boundary Extension (STARWALL, CARIDDI)"
grand_parent: "Model Extensions"
layout: default
render_with_liquid: false
---

## Stellarator no-wall boundary (model 183)

CARIDDI can provide the plasma boundary response for the stellarator model 183 without a conducting wall. CARIDDI needs the vacuum field $\nabla\chi$ on its control surface, which it obtains by calling [jorek2_gradchi_domm](/JOREK/diagnostics/jorek2_gradchi.html). Build that diagnostic with the same hard-coded parameters as the JOREK run.

The response file written by CARIDDI must be renamed to `starwall-response.dat`. It has to be computed for the same boundary grid and the same toroidal modes (`n_tor`, `n_period`) as the JOREK run.

Namelist settings for the JOREK run:

    CARIDDI_mode    = .true.    ! read the CARIDDI response format
    freeboundary    = .true.    ! use the boundary response
    resistive_wall  = .false.   ! no wall, ideal response only
    bc_natural_open = .false.
