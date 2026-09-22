---
title: "Free-boundary diagnostics and tools"
nav_order: 6
parent: "Free Boundary Extension (STARWALL, CARIDDI)"
grand_parent: "Model Extensions"
layout: default
render_with_liquid: false
---

This page is intended to describe tools related to the free-boundary extension. Here you will find

* [Construct a free-boundary equilibrium without knowing the coil currents](#construct-a-free-boundary-equilibrium-without-knowing-the-coil-currents)
* [Active controller model for vertical stabilization](#active-controller-model-for-vertical-stabilization)
* [Calculate wall forces](#calculate-wall-forces)
* [Calculate magnetic fields in vacuum](#calculate-magnetic-fields-in-vacuum)
* [Compress STARWALL response matrices](#compress-starwall-response-matrices)


## Construct a free-boundary equilibrium without knowing the coil currents

The previous section assumes that you already know the right coil geometries and currents required for your free boundary equilibrium. But what if you just want to create a new plasma in JOREK and you want to perform free-boundary studies (like disruptions or VDEs)? Also, what if you want to mock-up ferromagnetic structures (like the JET iron core) by varying coil currents?

We have implemented a routine that simplifies this process, eliminating the need to manually test different coil currents.

### Steps to follow

1. Create the fixed boundary equilibrium that you want to convert into free-boundary.
2. Copy the `jorek00000.h5` file to `jorek_restart.h5` (to use that plasma equilibrium).
3. In the JOREK code folder, run: `make -j 8 jorek2_find_Icoils`
4. Prepare a coil geometry file. Example files are available:
   - [aug_coils.zip](../assets/freebound/aug_coils.zip)
   - [iter_coils.zip](../assets/freebound/iter_coils.zip)
5. Name your coil geometry file `coils_geo.txt`.
6. Run the program: `./jorek2_find_Icoils < input_JOREK`
7. Your computed coil currents are stored in `Icoils_found.txt`.
8. **Important note**: The coil currents found are total currents (already multiplied by the number of turns). If you specified STARWALL coils with a number of turns, divide the found currents by the number of turns given by STARWALL before using them in the JOREK input file.

### For JET 

To predict the currents in the 10 JET circuits that feed the 20 JET coils:

1. Set `tokamak_device='JET'` in the JOREK input file.
2. Use the JET coils geometry file: [jet_coils.zip](../assets/freebound/jet_coils.zip).
3. To translate the found total JET currents to STARWALL format, divide each coil current by its number of turns. Use this [excel file](../assets/freebound/turns_cf_star.xlsx) to assist with the conversion.

## Active controller model for vertical stabilization
A PID controller was added to JOREK to enable active stabilization in free boundary simulation including n=0.

### Fundamentals

It acts on specific coils based on the vertical axis position and its evolution in the simulation according to:

$$\Delta I = K_P \, e(t) + K_D \frac{\mathrm{d} e(t)}{\mathrm{d} t} + K_I \int_{t_0}^{t} e(t)$$

where $e=Z(t)-Z_{\text{ref}}(t)$ is the deviation from the reference value and $K_P$, $K_I$, $K_D$ are the proportional, integral and derivative gains.

### Implementation in JOREK

The equation above is discretized and built into the coil_current_source routine in vacuum/vacuum_response.dat after the current profile of the PF coils is interpolated. There, it has the following form:

$$\begin{align}
\mathrm{d}Z_{\text{err}} &= Z_{\text{axis}}(t_n) - Z_{\text{ref}}(t_n)\\
\mathrm{d}Z_{\text{der}} &= \frac{Z_{\text{axis}}(t_n) - Z_{\text{axis}}(t_{n-1})}{\Delta t}\\
\mathrm{d}Z_{\text{integral}} &+= \left(Z_{\text{axis}}(t_n) - Z_{\text{ref}}(t_{n})\right)\Delta t \\
\Delta I &= K_P \, \mathrm{d}Z_{\text{err}} + K_D \, \mathrm{d}Z_{\text{der}} + K_I \, \mathrm{d}Z_{\text{integral}}\\
I(i) &= I(i) + \mathrm{vert\_FB\_amp\_ts}(i) \, \Delta I  
\end{align}$$

It has the following features with the parameters listed in the table below:
  * Activate the controller after a certain time start_VFB_ts
  * Use a constant reference value or an input profile
  * Specify tact time, a periodic interval after which the controller acts
  * Set coil current limits for each coil individually
  * Set an amplification factor for each coil individually to distribute the action

The controller is only activated if the amplification factors are specified and if the start time is exceeded.

### Parameters

| Parameter      | Comment   | Typical values           |
|---|---|---|
| vert_FB_gain(1)                 | proportional gain $K_P$                | 1.e0 (case dependent)  |
| vert_FB_gain(2)                 | **derivative** gain $K_D$              | 1.e0                   |
| vert_FB_gain(3)                 | **integral** gain $K_I$                | 1.e0                   |
| start_VFB_ts [t<sub>JOREK</sub>]| start time of vertical feedback               | -                      |
| vert_FB_amp_ts                  | gain amplification factor for individual coils| >0 for upper, <0 for lower coils    |
| I_coils_max [A]                 | maximum allowed current in the PF coil        | case specific     |
| vert_FB_tact [t<sub>JOREK</sub>]| Apply VFB only periodically                   | optional          |
| Z_ref_ts(t)                     | time trace of reference axis position         | case specific     | 
| vert_pos_file [t<sub>JOREK</sub> \| m] | input file for time dependent axis position   |                   |

### How to specify the input profile

The position can be specified by either 
  * an input profile file (vert_pos_file)
    * the time is given in [t<sub>JOREK</sub>]
    * the vertical position in [m]
    * when the simulation is longer than specified by the profile, the last value of the profile is used as the reference.  
  * the input parameter Z_axis_ref also used for the free boundary equilibrium
  * If none of this is specified, the initial or restart value of the axis is used

### Tuning

The optimal gains can vary from case to case depending on the equilibrium parameters and the use case.

Normally, the axis oscillates during the first time steps, when the controller is active. 
This oscillation can be used to adjust the gains. Otherwise, a sudden variation in the target value (a step response) can be used to improve the settings. 

The tuning can be performed in three steps:
  - Increase $K_P$ with the other gains 0 until the oscillation amplitude is low and the overshoot is small in case of a step response. If $K_P$ is too high the response will overshoot and the system becomes unstable.
  - Increase $K_D$ to reduce the oscillations
  - Increase $K_I$ if there is an error in steady state. (I haven't seen a big effect of this so far and the error was usually very low with only $K_P$ and $K_D$ gain)



## Calculate wall forces

This program reads JOREK restart files and calculates the total wall forces. **You need to run your simulations in full free-boundary mode (including the n=0 mode)**. Otherwise you cannot know the vacuum fields outside the plasma domain. 

It involves expensive volume integrals in the plasma, so run it in parallel with several MPI processes for large cases. Submit the same job scripts that you would use for JOREK, replacing `jorek_modelXXX` with `jorek2_wall_forces`. To compile it, type:

```bash
make -j 8 jorek2_wall_forces
```

To select JOREK restart files, create a file named `wall_forces.nml` with the following format:

```fortran
&wall_forces
   istart     = 100   ! start index of restart file
   iend       = 1000  ! end index of restart file
   delta_step = 40    ! step interval
   scale_fact = 1.01  ! scaling factor of wall (1.01 default). 
   n_phi_int  = 64    ! number of points for integration in toroidal direction (for B-field)
/
```
The forces are exported to `total_wall_forces_istart..iend..scale_fact.dat`.

#### Legacy format

For older versions, `wall_forces.nml` requires a single line:

```
0 1000 100
```
The result is stored in `total_wall_forces.dat`.

### Important note

To compute the total wall forces, an integral over a closed surface containing the STARWALL wall is performed ([see equation 7 in this paper](https://iopscience.iop.org/article/10.1088/1741-4326/aa8876/meta)). The routine assumes a specific ordering of the triangles describing this surface (a scaled value of the STARWALL wall). Depending on how this ordering is done in STARWALL, the triangle normals may flip sign and thus affect the sign of the computed forces. **TODO: Enforce always outward-pointing normals in the routine to ensure the correct sign of the total forces.**

## Calculate magnetic fields in vacuum 

This program reads JOREK restart files and calculates the magnetic field at arbitrary (x,y,z) points. The particularity with respect to `jorek2_postroc` is that you can get the fields also outside the JOREK grid. For example, it could be used to produce GEQDSK files from JOREK restart files.

To get the complete fields **you need to run your simulations in full free-boundary mode (including the n=0 mode)**. Otherwise you cannot know the vacuum fields outside the plasma domain. If you are running with `freeboundary=.false.` then you only get the fields produced by the plasma currents.

It involves the computation of expensive volume integrals in the plasma, so run it in parallel with several MPI processes for large cases. Submit the same job scripts that you would use for JOREK, replacing `jorek_modelXXX` with `jorek2_fields_xyz`. To compile it, type:

```bash
make -j 8 jorek2_fields_xyz
```

As input you must provide the points in cartesian coordinates in the file `xyz.nml` in the format:

```
number_of_points
x_1       y_1       z_1  
x_2       y_2       z_2  
 :         :         :   
 :         :         :   
x_n       y_n       z_n  
```

**DO NOT** specify the points exactly at the STARWALL coil/wall triangles. Otherwise singularities at those points may occur!

To select the given JOREK restart files create a file named `steps.nml` with 3 integers in the first line (`istart`, `iend`, `delta_step`). For example, to read the restart files from 100 up to 520 with intervals of 10 steps, write this first line in `steps.nml`:

```
100  520  10
```

The fields are exported in the file `fields_xyz.dat` and are also separated into coil, wall, and plasma contributions. 

**Note:** To calculate the required plasma fields, `n_plane` typically should be more than 60 (even for 2D).

## Compress STARWALL response matrices

> **Warning:** Response-matrix compression has not been merged into the standard `develop` branch. The current instructions apply only to the unmerged feature branch and may change before integration.

An experimental tool can compress STARWALL response matrices using a singular-value decomposition. See [Compress the Response Matrices](compress_response_matrices.html) for the feature-branch requirements, compilation instructions, parameters, and current limitations.
