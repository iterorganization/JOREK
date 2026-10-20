---
title: "Diagnostics and Scripts"
nav_order: 1
render_with_liquid: false
parent: "Diagnostics and post-processing"
has_children: true
nav_fold: true
---

This page lists some of the diagnostics available to JOREK users. This list is far from complete, and you can help by expanding it.

## Diagnostics

There are a lot of diagnostics available for JOREK. Some of the most important:

| Program | Description |
| :--- | :--- |
| **[jorek2vtk(_3d)](/JOREK/howto/introduction_to_jorek_diagnostics.html#jorek2vtk)** | Convert `jorek_restart.h5` files to VTK files |
| **[jorek\_read\_h5](jorek_read_h5.html)** | Python script for reading HDF5 files to VTK objects and files |
| **[jorek2vtk\_gaussvortterms](/JOREK/howto/introduction_to_jorek_diagnostics.html#beyond)** | Calculates the different terms of the vorticity equation integrated in the Gaussian points and creates a VTK file where the different terms can be visualized in the poloidal plane. |
| **[jorek2\_four](jorek2_four.html)** | Perform a two-dimensional Fourier analysis of the `jorek_restart.rst` file |
| **[jorek2\_poincare](/JOREK/howto/introduction_to_jorek_diagnostics.html#jorek2_poincare)** | Create a Poincaré plot for a JOREK restart file |
| **[jorek2\_postproc](/JOREK/howto/introduction_to_jorek_diagnostics.html#jorek2_postproc)** | The diagnostic program jorek2_postproc is a very flexible diagnostic tool that can be used interactively or using small scripts |
| **[jorek2\_diagno](/JOREK/howto/introduction_to_jorek_diagnostics.html#beyond)** | Extracts some useful data from a JOREK restart file |
| **jorek2\_diagno\_spi** | Extracts some useful data from a JOREK restart file, including the SPI fragments position and radius at the time of the restart file |
| **[jorek2\_connection2](/JOREK/howto/introduction_to_jorek_diagnostics.html#beyond)** | This diagnostic produces a poincaré plot of the magnetic field geometry |
| **[jorek2\_connection\_flux\_aligned](jorek2_connection_flux_aligned.html)** | Calculates connection length to the simulation boundary and the corresponding strike points of field lines on n=0 flux surfaces |
| **[jorek2\_fields\_xyz](jorek2_fields_xyz.html)** | Calculates B_x, B_y, B_z at given (x,y,z) cartesian coordinates |
| **[jorek2\_wall\_forces](jorek2_wall_forces.html)** | Calculates the total wall forces (F_x, F_y, F_z). It needs [JOREK-STARWALL](/JOREK/howto/running_jorek-starwall_for_the_first_time.html) |
| **[jorek2\_fast\_camera](jorek2_fast_camera.html)** | This diagnostic produces a camera image of the visible plasma (Dalpha) |
| **[rst\_bin2hdf5](rst_bin2hdf5.html)** & **[rst\_hdf52bin](rst_bin2hdf5.html#rst_hdf52bin)** | Convert `jorek_restart.rst` to `jorek_restart.h5` and vice-versa. |
| **[plot\_live\_data.sh](plot_live_data.sh.html)** | Plot several 0D quantities as the simulation runs (gnuplot) |
| **[Diagnostic Framework new\_diag](/JOREK/howto/introduction_to_jorek_diagnostics.html#new_diag)** | Allows to **evaluate arbitrary physical expressions ([current list of expressions](/JOREK/howto/introduction_to_jorek_diagnostics.html#expressions-available)) at arbitrary positions** in the computational domain either in **JOREK units or in SI units**. This framework is used by the [jorek2\_postproc](/JOREK/howto/introduction_to_jorek_diagnostics.html#jorek2_postproc) diagnostics for most of its functionality, but even more powerful diagnostics can easily be written on top of the framework. |
| **[radiation\_function\_diagno](/JOREK/physics/model_extensions/adas.html#checking-the-radiation-power-function-and-its-components)** | Output the radiation power function for the hard-coded impurity under hard-coded parameters |
| **[jorek2\_fieldlines\_vtk\_newdiag](jorek2_fieldlines_vtk_newdiag.html)** | This diagnostic evaluates arbitrary physical expressions (using the new_diag framework) along field lines and writes this information to both vtk and txt files |
| **[jorek2\_solcurrent](jorek2_solcurrent.html)** | This diagnostic postprocesses `jorek_restart.h5` to determine scrape-off layer currents (thermo-electric current between target plates) by solving [Stangeby](https://doi.org/10.1201/9780367801489) eq. 17.29 |
| **[fluid_loads_on_3D_wall](fluid_wall_load.html)** | Calculates fluid fluxes (heat and current loads) on a 3D thin walls discretized with triangles as explained in [ArtolaNF2026](https://iopscience.iop.org/article/10.1088/1741-4326/ae80e3)  |
| **[free-boundary diagnostics](/JOREK/physics/model_extensions/freebound/freebnd_tools.html)** | Here you can find diagnostics related to the free-boundary extension, such as the calculation of the B-field outside the grid and the computation of wall forces |

The diagnostics tools should have been compiled with the same hard-coded parameters as JOREK.

To **run diagnostics on Marconi**, you may need to book a debug node with a command of the type:

```
salloc --nodes=1 --time=01:00:00 --exclusive --account=FUA36_MHD --partition=skl_fua_dbg
```

To **run diagnostics on Leonardo DCGP**, you may need to book a debug node with a command of the type:

```
salloc --nodes=1 --time=10:00 --account=FUA38_MHD_0 --partition=dcgp_fua_dbg
```

and run:

```
srun jorek2_diagno < jorek_in
```

Alternatively you can run directly on the login node by first setting the ulimit memory limit to unlimited with:

```
ulimit -s unlimited
```

## Scripts

A couple of scripts are available in the folder **trunk/util/**:

| Script | Description |
| :--- | :--- |
| **[config.sh](hard-coded_parameters.html)** | List or modify **hard-coded parameters** |
| **[plot\_live\_data.sh](plot_live_data.sh.html)** | **Plot time-traces of energies or growth rates** (and more) while a simulation is running or afterwards |
| **[plot\_grids.sh](plot_grids.sh.html)** | **Plot the grid to a file or the screen** |
| **[plot\_live\_data.py](plot_live_data.py.html)** | **Python version of plot_live_data.sh** |
| **[read\_jorek\_logfile.py](read_jorek_logfile.py.html)** | **Search JOREK logfile for specific variable** |
| **[plot\_mlog.py](plot_mlog.py.html)** | **Plot memory usage** |

## Output files

During a run some output files are generated. This is a short summary of the different files.

- [macroscopic\_vars](macroscopic_vars.html) contains a list of the energies and growth rates in the simulation (and various other quantities)
- `jorek00XXX.rst` are snapshots that can be used to restart the simulation
- `boundary.txt`
- `equilibrium.txt`
- `grid_fluxsurface.dat`
- `grid_initial.dat`
- `jorek2.ps` is generated at the end of the simulation and contains some plots of the energies and magnetic field.
- `qprofile.dat`
- `special_equilibrium_points.dat`
- `T_rho_profiles.dat`
