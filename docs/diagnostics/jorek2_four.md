---
title: "jorek2_four"
nav_order: 28
render_with_liquid: false
parent: "General"
---

## jorek2_four

The program `jorek2_four` (`make jorek2_four`) allows to perform a
two-dimensional Fourier analysis of the data contained in the file
`jorek_restart.rst`. The Fourier analysis is performed in magnetic
coordinates.

      jorek2_four < jorek-inputfile

Some numerical parameters can be modified by putting a namelist file `four_params.nml` in the current directory. The namelist file with the most important parameter and default values:

              &four_params
                nstpts      = 30   ! Number of radial points
                deltaphi    = 0.3  ! Step size in toroidal direction
                rad_range = 0.001, 0.999 ! Radial range of normalized psi
              /

The parameter `deltaphi` should be chosen such that about 200 to 1000 points are obtained per flux surface from field line tracing. There are some additional parameters which usually don’t have to be modified (refer to the code for details).

- The output is written to ascii files `Flux_modes_n000` and can easily be plotted using `gnuplot` or `xmgrace`.

- The main program is implemented in `diagnostics/jorek2_four.f90` and the core routines are found in `diagnostics/mod_fourier.f90`.

- Implementation details:
  - The code is OpenMP parallelized.
  - First, field line tracing in the axisymmetric component of the magnetic field is performed with the number of start points given by the parameter `nstpts`.
  - A fixed number of poloidal points (determined from the poloidal mode numbers `m_pol_range`) is determined for each flux surface. These points are equidistant in the poloidal magnetic angle $\theta_{mag}$.
  - The physical variables are Fourier-transformed after being interpolated to these equidistant points. FFTW is used for the transformation.
  - With `debug = .true.` in the namelist, additional debugging output is created which helps to search for problems.
