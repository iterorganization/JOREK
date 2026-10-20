---
title: "JOREK-STARWALL tutorial: resistive-wall tearing mode"
nav_order: 2
parent: "Free Boundary Extension (STARWALL, CARIDDI)"
grand_parent: "Model Extensions"
layout: default
render_with_liquid: false
---

# Running JOREK-STARWALL for the first time

This tutorial shows how to run a small JOREK simulation coupled to STARWALL. It uses the `model600/intear_freebnd` tearing-mode case and replaces the ideal conducting boundary condition with the response of a simple resistive wall.

Before starting, you should be able to compile and run JOREK as described in [Running JOREK for the first time](/JOREK/howto/running_jorek_for_the_first_time.html). You also need a compiled STARWALL executable; follow [Getting the STARWALL source code and compiling it](running_STARWALL.html#getting-the-source-code) rather than using a machine-specific configuration copied from this tutorial.

The coupling is prepared in three stages:

1. Run JOREK once with `freeboundary=.t.` to write `boundary.txt`.
2. Run STARWALL with that boundary to generate `starwall-response.dat`.
3. Run the JOREK time evolution with `starwall-response.dat` in the run directory.

STARWALL is needed only when generating the response file, not while JOREK advances the time evolution.

## 1. Prepare the JOREK case

Choose paths for the JOREK source tree and the run directory. The example below keeps the run outside the source tree:

```bash
JOREK_DIR=/absolute/path/to/JOREK
RUN_DIR=/absolute/path/to/intear2_reswall

mkdir -p "$RUN_DIR"
cp "$JOREK_DIR/namelist/model600/intear_freebnd" "$RUN_DIR/intear"

cd "$JOREK_DIR"
# choose minimalistic model options for faster runs (no impurities, no vpar, only one T, constant density)
./util/config.sh model=600 n_tor=3 n_period=1 n_plane=4 with_impurities=.false. with_vpar=.false. with_TiTe=.false. with_rho=.false.

make clean
make -j 8
make -j 8 jorek2vtk jorek2_postproc

cp jorek_model600 jorek2vtk jorek2_postproc "$RUN_DIR/"
ln -sfn "$JOREK_DIR/util" "$RUN_DIR/util"
```

For this configuration, JOREK contains the axisymmetric component and the cosine and sine components of the $n=1$ harmonic. The STARWALL input below supplies a free-boundary response only for $n=1$; the $n=0$ component therefore remains fixed boundary.

All JOREK executables used with the same restart files must be compiled with consistent hard-coded parameters.

## 2. Copy the STARWALL executable

After following the [STARWALL compilation instructions](running_STARWALL.html#compiling), copy its executable into the run directory:

```bash
STARWALL_DIR=/absolute/path/to/STARWALL
cp "$STARWALL_DIR/o_3d/STARWALL_JOREK_Linux" "$RUN_DIR/"
```

If your STARWALL configuration places the executable in a system-specific subdirectory, adjust the source path accordingly.

## 3. Generate `boundary.txt` with JOREK

Set `freeboundary = .t.` in the input file (intear), then, run JOREK once:

```bash
export OMP_NUM_THREADS=8
mpirun -n 2 ./jorek_model600 < intear | tee logfile_boundary
```

This first invocation is expected to stop because `starwall-response.dat` does not exist yet. Before stopping, JOREK should write the plasma coupling surface to `boundary.txt`. Confirm that the file is present and non-empty:

It is important to generate this file with `freeboundary=.t.`. A boundary exported by a fixed-boundary calculation can contain a different number of boundary elements and must not be reused for this step. The general response-file workflow is summarized in [Getting the response file](jorek_freebnd_params.html#getting-the-response-file).

 **Note:** JOREK provides several input parameters for configuring the free-boundary extension, including the calculation of free-boundary equilibria (the $n=0$ mode). See [Free-boundary input parameters](jorek_freebnd_params.html) for details.

## 4. Create the STARWALL input

Create a file named `input_starwall` in the run directory:

```fortran
&PARAMS
  ! Calculate the response for a resistive wall.
  ! 0 = no wall, 1 = ideal wall, 2 = resistive wall
  i_response = 2,

  ! Include only the n=1 toroidal harmonic.
  n_harm = 1,
  n_tor  = 1,

  ! Toroidal grid points in the real-space representation.
  nv = 40,

  ! Numerical parameters normally left unchanged.
  delta    = 0.001,
  n_points = 10,
  nwall    = 1,

  ! Represent the wall with a Fourier series.
  iwall = 1,
/

&PARAMS_WALL
  ! Thin-wall resistivity in ohms.
  eta_thin_w = 1.d-4,

  ! Poloidal and toroidal wall-grid points.
  nwu = 32,
  nwv = 32,

  ! R = 10 + 1.10 cos(theta), Z = 1.10 sin(theta).
  mn_w = 2,
  n_w  =  0,    0,
  m_w  =  0,    1,
  rc_w = 10., 1.10,
  rs_w =  0., 0.00,
  zc_w =  0., 0.00,
  zs_w =  0., 1.10,
/
```

This input describes a circular wall centred at $R=10$, $Z=0$, with minor radius $1.10$. For the meaning of the parameters, wall-resistivity units, alternative wall representations, and coil inputs, consult the maintained [STARWALL parameter setup](running_STARWALL.html#parameter-setup).

The names `n_harm` and `n_tor` can be confusing: in STARWALL, `n_harm` is the number of toroidal harmonics in the response, while `n_tor` lists their mode numbers. These are not the same as the JOREK compile-time `n_tor`. See [Select modes and harmonics](running_STARWALL.html#select-modes-and-harmonics) for examples.

## 5. Generate `starwall-response.dat`

Run STARWALL from the directory containing both `boundary.txt` and `input_starwall`:

```bash
cd "$RUN_DIR"
export OMP_NUM_THREADS=1
mpirun -n 16 ./STARWALL_JOREK_Linux ./input_starwall | tee logfile_starwall
```

STARWALL is MPI-only, so use one OpenMP thread per MPI task. The calculation should take a few minutes for this small case. For a batch system, adapt the resource request to the machine as described under [Submitting the STARWALL job](running_STARWALL.html#submitting-the-job); old machine-specific job scripts should not be copied unchanged.

STARWALL also writes `wall.vtk`, which contains the triangular discretization of the wall and can be inspected with ParaView or VisIt.

## 6. Run the coupled JOREK time evolution

The run directory now contains the required response file, and `freeboundary=.t.` is already set in `intear`. Start JOREK in the same way as a standalone calculation:

```bash
cd "$RUN_DIR"
export OMP_NUM_THREADS=8
mpirun -n 2 ./jorek_model600 < intear | tee logfile_jorek
```

For a production run, use the batch submission procedure appropriate to your HPC system.

## 7. Check and analyze the run

Near the beginning of `logfile_jorek`, the `STARWALL RESPONSE INFORMATION` block should report the response dimensions, wall resistivity, and included harmonics. For this example, the harmonic list should include the cosine and sine components of $n=1$:

```text
i_tor = 2 (n=1 cos), 3 (n=1 sin)
```

JOREK writes `wallcurr.XXXXXX.vtk` at the interval controlled by `nout`. These files contain the wall-current stream function and can be opened with ParaView or VisIt. General diagnostic and plotting instructions are available in [Introduction to JOREK diagnostics](/JOREK/howto/introduction_to_jorek_diagnostics.html).

To study the physical effect of the wall:

1. Run otherwise identical fixed- and free-boundary cases and compare their linear growth rates and mode structures.
2. Vary the JOREK input parameter `wall_resistivity_fact` to scale the wall resistivity without regenerating the STARWALL response. This can be used to explore the limits from an ideal wall to a highly resistive wall. See the [free-boundary input parameters](jorek_freebnd_params.html#main-parameters).
3. Repeat the calculation with $n=1$ and $n=2$. Reconfigure JOREK so that both harmonics and sufficient toroidal planes are present, then set `n_harm=2` and `n_tor=1,2` in `input_starwall`. The JOREK and STARWALL harmonic sets must be consistent.

## Optional: add a diagnostic coil

Diagnostic coils are defined when generating the STARWALL response. Add the following entry to the `&PARAMS` block in `input_starwall`:

```fortran
diag_coil_file = 'diag_coils.nml',
```

Then create `diag_coils.nml` in the run directory:

```fortran
&coil_set_nml
  description = 'example coil set'
  ncoil       = 1
/

&coils_nml
  coil(1)%name      = 'diag01'
  coil(1)%coil_type = 'general_thin'
  coil(1)%resist    = 1.d0
  coil(1)%n_pts     = 5
  coil(1)%width     = 0.01
  coil(1)%xpts      = 11.,  11., 11.,  11., 11.
  coil(1)%ypts      = -0.5, 0.5, 0.5, -0.5, -0.5
  coil(1)%zpts      = -0.5,-0.5, 0.5,  0.5, -0.5
  coil(1)%nturns    = 1
/
```

Run STARWALL again to regenerate `starwall-response.dat`, then rerun JOREK. Plot the diagnostic-coil signal with:

```bash
./util/plot_live_data.sh -q diag_coil_curr
```

The supported coil types and all coil fields are documented under [Include coils](running_STARWALL.html#include-coils).

## Next steps

For more advanced analysis, see the [free-boundary diagnostics and tools](freebnd_tools.html), including wall-force and vacuum-field calculations. The papers that should be cited for JOREK-STARWALL work are listed under [Free-boundary references](references.html).
