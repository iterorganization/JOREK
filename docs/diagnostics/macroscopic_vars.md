---
title: "Macroscopic_vars.dat"
nav_order: 25
render_with_liquid: false
parent: "Diagnostics and post-processing"
---

This file contains the energies and growth rates at specific time steps in the simulation, and can be used to plot the variation of these quantities in time with the [plot_live_data.sh](plot_live_data.sh.md) script.

The format is 
   @rcs_version:  GIT revision        b575260-dirty
   @jorek_model:   303
   @n_tor:     3
   @n_plane:     8
   @n_period:     1

And then a series of values at specific times
   @times:     1   1.000000000E+03
   @energies:  1.000000000E+03  2.012235987E-02  2.073465471E-22  1.291279331E-16  1.245421010E-19
   @growth_rates:  1.500000000E+03 -8.745715263E-13  6.921484349E-04  7.638978557E-06 -2.518593016E-06
