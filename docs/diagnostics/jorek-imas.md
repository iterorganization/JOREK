---
title: "JOREK-IMAS"
nav_order: 5
parent: "Diagnostics and post-processing"
layout: default
render_with_liquid: false
---


# JOREK-IMAS

## What Is IMAS?

The Integrated Modelling & Analysis Suite (IMAS) provides a common data model and a software ecosystem for exchanging experimental and simulation data between fusion codes. The [IMAS Data Dictionary](https://github.com/iterorganization/IMAS-Data-Dictionary) implements this model as machine-independent **Interface Data Structures (IDSs)** with standardized names and structures.

**New users should start with the [open-source IMAS tutorial](https://github.com/iterorganization/IMAS-tutorial).** Its notebooks introduce IMAS, the Data Dictionary, data access, IDS exploration and visualization, and simulation-database management. The tutorial is a living document and links to the maintained documentation for each component.

## Open-Source IMAS Components

Most of the IMAS software relevant to JOREK users is developed openly:

| Component | Purpose |
| --- | --- |
| [IMAS Data Dictionary](https://github.com/iterorganization/IMAS-Data-Dictionary) | Defines the IMAS data model and IDS schemas for experimental and simulation data. |
| [IMAS-Python](https://github.com/iterorganization/IMAS-Python) | Python library for creating, reading, writing, and manipulating IDSs. |
| [IMAS-Fortran](https://github.com/iterorganization/IMAS-Fortran) | Fortran data-access library used to manipulate IDS structures and connect them to the available storage backends. |
| [IMAS-ParaView](https://github.com/iterorganization/IMAS-ParaView) | ParaView plugins for visualizing GGD meshes and other IMAS data. This replaces the former GGD-VTK workflow. |
| [IMAS-MUSCLE3](https://github.com/iterorganization/IMAS-MUSCLE3) | Integration tools for coupling IMAS-compatible codes in MUSCLE3 workflows. Its actors exchange IDS data between the coupled components. |
| [SimDB](https://github.com/iterorganization/SimDB) | Command-line and client-server tools for tracking, querying, transferring, and managing simulations and their metadata. |
| [SimDB Dashboard](https://github.com/iterorganization/SimDB-Dashboard) | Web interface for querying and viewing simulations in a remote SimDB database. |

Many more open-source fusion and IMAS repositories are available from the [ITER Organization on GitHub](https://github.com/iterorganization).

## How Can JOREK Users Benefit from IMAS?

- **Obtain initial equilibria:** Build a JOREK input file from an equilibrium IDS, including $FF'$, coil currents, and the poloidal-flux map. Use [SimDB](https://github.com/iterorganization/SimDB) or the [SimDB Dashboard](https://github.com/iterorganization/SimDB-Dashboard) to find available simulations and their database metadata.

- **Obtain tokamak components:** Read specifications such as vacuum-vessel geometry, first-wall panels, diagnostics, and coils from the corresponding IDSs when they are included in an accessible IMAS entry.

- **Archive, share, and reference simulations:** Export JOREK results to IDSs and manage the resulting simulation and its metadata with SimDB.

- **Visualize and process IMAS data:** Use [IMAS-ParaView](https://github.com/iterorganization/IMAS-ParaView) for GGD meshes, machine descriptions, and diagnostic geometries, or work directly with IDSs through [IMAS-Python](https://github.com/iterorganization/IMAS-Python) and [IMAS-Fortran](https://github.com/iterorganization/IMAS-Fortran).

- **Couple JOREK to other codes and create workflows:** MUSCLE3 can couple JOREK to other codes, for example to calculate synthetic diagnostic signals from JOREK output. [IMAS-MUSCLE3](https://github.com/iterorganization/IMAS-MUSCLE3) provides the actors that exchange IDS data between the coupled components.

## Creating a JOREK Input File from an IMAS Entry

1. Follow the [IMAS tutorial](https://github.com/iterorganization/IMAS-tutorial) to become familiar with IDSs and the current data-access tools.
2. Find an appropriate simulation with the [SimDB Dashboard](https://github.com/iterorganization/SimDB-Dashboard). Record the user, database, pulse, and run identifying the IMAS entry, then select a time slice containing a suitable equilibrium.
3. Configure the IMAS data-access environment required at your site. Installation and usage information is maintained in the [IMAS-Python](https://github.com/iterorganization/IMAS-Python) repository; module names and available storage backends may be site-specific.
4. From `communication/IMAS`, run `imas2jorek.py`, replacing the uppercase placeholders with the selected entry and time:

```bash
python imas2jorek.py \
    --user USER \
    --database DATABASE \
    --pulse PULSE \
    --run RUN \
    --time TIME
```

Display all arguments and defaults, including the storage backend, Data Dictionary major version, IDS occurrence, and tokamak used to construct the boundary, with:

```bash
python imas2jorek.py -h
```

The script operates similarly to `eqdsk2jorek`. The `--tokamak` (`-tk`) option selects the tokamak used to construct the initial $R$-$Z$ boundary; inspect and adapt the generated grid boundary as needed.

After `imas2jorek.py` finishes, it produces:

- `jorek_namelist`, containing the poloidal-flux boundary conditions and coil currents.
- `jorek_density`, `jorek_temperature`, and `jorek_ffprime`, containing the profiles.

## Exporting JOREK Data to IMAS

Load the IMAS module. On the ITER cluster:

```bash
module load IMAS
```

The `jorek2_IDS` tool converts JOREK HDF5 restart files into IMAS data. Enable the IMAS interface and provide the Access Layer compiler and linker flags in `Makefile.inc`:

```makefile
USE_IMAS = 1

IMASINCLUDE=`pkg-config --cflags-only-I al-fortran`
IMASLIB=`pkg-config --libs-only-L --libs-only-l al-fortran`

# For the old Access Layer 4, use these lines instead:
# IMASINCLUDE=`pkg-config --cflags-only-I imas-ifort`
# IMASLIB=`pkg-config --libs-only-L --libs-only-l imas-ifort`
```

Configure `jorek2_IDS` consistently with the model and toroidal resolution used by the simulation, then compile it. For example:

```bash
util/config.sh n_tor=1 n_period=1 n_plane=1 with_TiTe=.false. with_vpar=.false. with_impurities=.true.
make -j 8 jorek2_IDS
```

In the simulation directory containing the restart files, create `imas.nml`. The comments below describe its parameters:

```fortran
&imas_params

  ! --- Where is this simulation going to be stored?
  user        = 'artolaj'           ! your username
  database    = 'test_database'     ! your local database, if it does not exist it will create a new one
  shot_number = 111111              ! choose a shot or reference number
  run_number  = 1                   ! choose a run number
  overwrite_entry = .false.         ! If true:  it will overwrite pre-existing entries in the database (CAREFUL!!)
                                    ! If false: it will append the data to existing IDSs

  ! --- Which restart files are going to be exported?
  i_begin     = 700           ! Do not consider restart files before 700
  i_end       = 5000          ! Do not consider restart files after 5000
  i_jump_steps= 100           ! jump restart files by this step

  ! --- Which IDSs are going to be exported?
  export_JOREK_variables  = .true.     ! --- Exports JOREK raw variables to occurrence 1 of plasma_profiles
  export_1d_profiles      = .false.    ! --- Exports 1D flux-averaged profiles to occurrence 0 of plasma_profiles
  export_equilibrium      = .false.
  export_summary          = .false.
  export_disruption       = .false.
  export_wall             = .false.    ! --- Needs free-boundary mode
  export_pf_passive       = .false.    ! --- Needs free-boundary mode
  export_pf_active        = .false.    ! --- Needs free-boundary mode
  export_radiation        = .false.    ! --- Only working with SPI simulations with marker model
  export_spi              = .false.    ! --- Needs SPI simulation
  export_field_extension  = .false.    ! --- Export psi(n=0) and B(3D) to a rectangular grid that goes beyond
                                       ! ---    the JOREK grid. The n=0 components are stored in
                                       ! ---    equilibrium/profiles_2D and the 3D fields in plasma_profiles
                                       ! ---    occurrence 2.

  ! --- Geometry of coils and passive components (as in STARWALL input)
  ! --- These are needed to export the coil geometries in the IDSs pf_active & pf_passive
  passive_coil_geo_file = 'passive_components.nml'
  active_coil_geo_file  = 'polcoils_iter.nml'

  ! --- If true, read only radiation projections for the non-coronal equilibrium model
  rad_only_projections_h5 = .false.

  ! --- Brief description of the simulation
  simulation_description = 'JOREK 2D disruption simulation.'

  ! --- Parameters defining the rectangular grid used in equilibrium and plasma_profiles occurrence 2.
  rect_grid_params%nR    = 100
  rect_grid_params%nZ    = 200
  rect_grid_params%R_min = 3.d0
  rect_grid_params%R_max = 10.d0
  rect_grid_params%Z_min = -6.d0
  rect_grid_params%Z_max = 6.d0

/
```

Run the converter with one MPI process:

```bash
mpirun -n 1 ./jorek2_IDS < jorek_input_file
```

Successful exports print messages such as `XXXXX IDS exported`. The stored data can then be inspected with [IMAS-ParaView](https://github.com/iterorganization/IMAS-ParaView) or processed directly with [IMAS-Python](https://github.com/iterorganization/IMAS-Python) and [IMAS-Fortran](https://github.com/iterorganization/IMAS-Fortran).

## Available IDS Nodes

The IDSs and variables that can currently be filled are listed below.

### `summary`

```text
summary_ids.ids_properties.homogeneous_time
summary_ids.ids_properties.comment
summary_ids.time
summary_ids.global_quantities.r0.value
summary_ids.global_quantities.b0.value
summary_ids.global_quantities.beta_pol.value
summary_ids.global_quantities.beta_tor.value
summary_ids.global_quantities.beta_tor_norm.value
summary_ids.global_quantities.ip.value
summary_ids.global_quantities.li_3.value
summary_ids.global_quantities.volume.value
summary_ids.global_quantities.q_95.value
summary_ids.global_quantities.energy_thermal.value
summary_ids.global_quantities.energy_b_field_pol.value
summary_ids.global_quantities.power_ohm.value
summary_ids.global_quantities.power_radiated.value
summary_ids.heating_current_drive.power_additional.value
summary_ids.local.magnetic_axis.position.r
summary_ids.local.magnetic_axis.position.z
summary_ids.local.magnetic_axis.position.psi
summary_ids.boundary.type.value
summary_ids.boundary.minor_radius.value
summary_ids.boundary.elongation.value
summary_ids.boundary.triangularity_upper.value
summary_ids.boundary.triangularity_lower.value
summary_ids.boundary.geometric_axis_r.value
summary_ids.boundary.geometric_axis_z.value
```

### `plasma_profiles` (occurrence 0)

Provides poloidal-flux-averaged profiles using the $n=0$ poloidal flux.

```text
plasma_profiles_ids.profiles_1d.grid.psi_magnetic_axis
plasma_profiles_ids.profiles_1d.grid.psi_boundary
plasma_profiles_ids.profiles_1d.grid.rho_pol_norm
plasma_profiles_ids.profiles_1d.grid.psi

plasma_profiles_ids.profiles_1d.t_i_average
plasma_profiles_ids.profiles_1d.electrons.temperature
plasma_profiles_ids.profiles_1d.electrons.density
plasma_profiles_ids.profiles_1d.pressure_thermal
plasma_profiles_ids.profiles_1d.phi_potential
plasma_profiles_ids.profiles_1d.conductivity_parallel
plasma_profiles_ids.profiles_1d.j_total
plasma_profiles_ids.profiles_1d.e_field.parallel
plasma_profiles_ids.profiles_1d.e_field.radial
plasma_profiles_ids.profiles_1d.ion.velocity.parallel
plasma_profiles_ids.profiles_1d.ion.velocity.poloidal
plasma_profiles_ids.profiles_1d.ion.velocity.diamagnetic
plasma_profiles_ids.profiles_1d.zeff
plasma_profiles_ids.profiles_1d.ion.density
plasma_profiles_ids.profiles_1d.ion.element.a
plasma_profiles_ids.profiles_1d.neutral.density
plasma_profiles_ids.profiles_1d.neutral.element.a
plasma_profiles_ids.profiles_1d.ion.element.z_n
plasma_profiles_ids.profiles_1d.q
plasma_profiles_ids.profiles_1d.grid.rho_tor_norm
```

### `plasma_profiles` (occurrence 1)

Contains high-precision raw JOREK variables. This IDS can be used to calculate current density and magnetic field accurately. It also contains electron and ion temperatures, mass density, electrostatic potential, and toroidal vorticity. Finite-element and Fourier-series coefficients are stored in GGD format and can be visualized with [IMAS-ParaView](https://github.com/iterorganization/IMAS-ParaView).

```text
plasma_profiles_ids.ids_properties.homogeneous_time
plasma_profiles_ids.time

plasma_profiles_ids.ggd.time
plasma_profiles_ids.ggd.psi
plasma_profiles_ids.ggd.phi_potential
plasma_profiles_ids.ggd.r_j_total_phi
plasma_profiles_ids.ggd.vorticity_over_r
plasma_profiles_ids.ggd.mass_density
plasma_profiles_ids.ggd.electrons.temperature
plasma_profiles_ids.ggd.t_i_average
plasma_profiles_ids.ggd.velocity_parallel_over_b_field
plasma_profiles_ids.ggd.ion.element.a
plasma_profiles_ids.ggd.ion.element.z_n
plasma_profiles_ids.ggd.ion.density
plasma_profiles_ids.ggd.neutral.element.a
plasma_profiles_ids.ggd.neutral.element.z_n
plasma_profiles_ids.ggd.neutral.density

##Grid
plasma_profiles.grid_ggd[0].identifier.description
plasma_profiles.grid_ggd[0].identifier.name
plasma_profiles.grid_ggd[0].identifier.index

plasma_profiles.grid_ggd[0].grid_subset[0].identifier.name
plasma_profiles.grid_ggd[0].grid_subset[0].identifier.index
plasma_profiles.grid_ggd[0].grid_subset[0].identifier.description
plasma_profiles.grid_ggd[0].grid_subset[0].dimension

plasma_profiles.grid_ggd[0].space[0].geometry_type.index
plasma_profiles.grid_ggd[0].space[0].coordinates_type
plasma_profiles.grid_ggd[0].space[0].objects_per_dimension[0].object[:].geometry
plasma_profiles.grid_ggd[0].space[0].objects_per_dimension[0].object[:].geometry_2d
plasma_profiles.grid_ggd[0].space[0].objects_per_dimension[2].object[:].nodes
plasma_profiles.grid_ggd[0].space[0].objects_per_dimension[2].object[:].geometry_2d

plasma_profiles.grid_ggd[0].space[1].coordinates_type
plasma_profiles.grid_ggd[0].space[1].geometry_type.index
plasma_profiles.grid_ggd[0].space[1].identifier.description
plasma_profiles.grid_ggd[0].space[1].objects_per_dimension[0].object[:].geometry
```

### `plasma_profiles` (occurrence 2)

Contains the $B_R$ and $B_Z$ magnetic-field components on a rectangular grid extending beyond the JOREK grid. It is available only for full free-boundary simulations with `export_field_extension = .true.`. The $(R,Z)$ space is a regular two-dimensional grid containing values at its nodes rather than Bézier coefficients, while the toroidal direction retains its Fourier representation. This occurrence is intended for field-line-tracing codes such as SMITER.

```text
plasma_profiles_ids.ids_properties.homogeneous_time
plasma_profiles_ids.time
plasma_profiles_ids.ggd.b_field

##Grid
plasma_profiles.grid_ggd[0].identifier.description
plasma_profiles.grid_ggd[0].identifier.name
plasma_profiles.grid_ggd[0].identifier.index

plasma_profiles.grid_ggd[0].grid_subset[0].identifier.name
plasma_profiles.grid_ggd[0].grid_subset[0].identifier.index
plasma_profiles.grid_ggd[0].grid_subset[0].identifier.description
plasma_profiles.grid_ggd[0].grid_subset[0].dimension

plasma_profiles.grid_ggd[0].space[0].geometry_type.index
plasma_profiles.grid_ggd[0].space[0].coordinates_type
plasma_profiles.grid_ggd[0].space[0].objects_per_dimension[0].object[:].geometry
plasma_profiles.grid_ggd[0].space[0].objects_per_dimension[2].object[:].nodes

plasma_profiles.grid_ggd[0].space[1].coordinates_type
plasma_profiles.grid_ggd[0].space[1].geometry_type.index
plasma_profiles.grid_ggd[0].space[1].identifier.description
plasma_profiles.grid_ggd[0].space[1].objects_per_dimension[0].object[:].geometry
```

### `wall`

The `wall` IDS contains two wall descriptions:

- The inner ITER vacuum-vessel layer, represented as a thin wall discretized with linear triangular elements. Current-density vector, effective thickness, and resistivity are provided for each triangle.
- An axisymmetric representation of the ITER first wall. Incident perpendicular heat flux, current-density vector for halo-current calculations, and poloidal flux are stored at the nodes of a triangular mesh.

These quantities can be visualized with [IMAS-ParaView](https://github.com/iterorganization/IMAS-ParaView).

```text
wall_ids.ids_properties.homogeneous_time
wall_ids.time

# Vacuum vessel inner layer
wall_ids.description_ggd[0].ggd[:].j_total
wall_ids.description_ggd[0].ggd[:].resistivity

# First wall surface
wall_ids.description_ggd[1].ggd[:].j_total
wall_ids.description_ggd[1].ggd[:].power_density
wall_ids.description_ggd[1].ggd[:].psi

# Grid
wall_ids.description_ggd[:].type.index
wall_ids.description_ggd[:].grid_ggd[0]identifier.index
wall_ids.description_ggd[:].grid_ggd[0]identifier.description
wall_ids.description_ggd[:].grid_ggd[0]grid_subset[0].identifier.index
wall_ids.description_ggd[:].grid_ggd[0]grid_subset[0].dimension
wall_ids.description_ggd[:].grid_ggd[0]space[0].identifier.index
wall_ids.description_ggd[:].grid_ggd[0]space[0].geometry_type.index
wall_ids.description_ggd[:].grid_ggd[0]space[0].coordinates_type
wall_ids.description_ggd[:].grid_ggd[0]space[0].objects_per_dimension[0].geometry_content.index
wall_ids.description_ggd[:].grid_ggd[0]space[0].objects_per_dimension[0].object[:].geometry
wall_ids.description_ggd[:].grid_ggd[0]space[0].objects_per_dimension[2].object[:].nodes

wall_ids.description_ggd[:].thickness[0].grid_subset[0].grid_index
wall_ids.description_ggd[:].thickness[0].grid_subset[0].grid_subset_index
wall_ids.description_ggd[:].thickness[0].grid_subset[0].values

wall_ids.description_ggd[:].component[0].type[0].grid_index
wall_ids.description_ggd[:].component[0].type[0].grid_subset_index
wall_ids.description_ggd[:].component[0].type[0].identifier.index
wall_ids.description_ggd[:].component[0].type[0].identifier.name
wall_ids.description_ggd[:].component[0].type[0].identifier.description
```

### `pf_passive`

The `pf_passive` IDS contains additional passive conductors in which currents are induced: the outer vacuum-vessel layer, outer triangular support, and divertor inboard rail. They are discretized as toroidal loops.

```text
pf_passive.ids_properties.homogeneous_time
pf_passive.time
pf_passive.loop.name
pf_passive.loop.resistance
pf_passive.loop.element.turns_with_sign
pf_passive.loop.element.geometry.geometry_type
pf_passive.loop.element.geometry.rectangle.r
pf_passive.loop.element.geometry.rectangle.z
pf_passive.loop.element.geometry.rectangle.width
pf_passive.loop.element.geometry.rectangle.height
pf_passive.loop.current
```

### `pf_active`

```text
pf_active.ids_properties.homogeneous_time
pf_active.time
pf_active.coil.name
pf_active.coil.resistance
pf_active.coil.element.turns_with_sign
pf_active.coil.element.geometry.geometry_type
pf_active.coil.element.geometry.rectangle.r
pf_active.coil.element.geometry.rectangle.z
pf_active.coil.element.geometry.rectangle.width
pf_active.coil.element.geometry.rectangle.height
pf_active.coil.current.data
```

### `disruption`

```text
disruption_ids.ids_properties.homogeneous_time
disruption_ids.time
disruption_ids.global_quantities.current_halo_pol
disruption_ids.global_quantities.current_halo_phi
disruption_ids.global_quantities.power_ohm
disruption_ids.global_quantities.power_ohm_halo
disruption_ids.global_quantities.power_parallel_halo
disruption_ids.global_quantities.power_radiated_electrons_impurities
```

### `spi`

```text
spi.ids_properties.homogeneous_time
spi.time
spi.injector.pellet.core.atoms_n
spi.injector.pellet.core.species.a
spi.injector.pellet.core.species.z_n
spi.injector.pellet.core.species.density
spi.injector.velocity_mass_centre_fragments_r
spi.injector.velocity_mass_centre_fragments_z
spi.injector.velocity_mass_centre_fragments_tor
spi.injector.fragment.position.r
spi.injector.fragment.position.z
spi.injector.fragment.position.phi
spi.injector.fragment.velocity_r
spi.injector.fragment.velocity_z
spi.injector.fragment.velocity_tor
spi.injector.fragment.volume
```

### `equilibrium`

```text
equilibrium_ids.ids_properties.homogeneous_time
equilibrium_ids.time

equilibrium_ids.time_slice.profiles_1d.psi
equilibrium_ids.time_slice.profiles_1d.pressure
equilibrium_ids.time_slice.profiles_1d.dpressure_dpsi
equilibrium_ids.time_slice.profiles_1d.f_df_dpsi
equilibrium_ids.time_slice.profiles_1d.j_parallel
equilibrium_ids.time_slice.profiles_1d.q
equilibrium_ids.time_slice.profiles_1d.rho_tor_norm

equilibrium_ids.vacuum_toroidal_field.r0
equilibrium_ids.vacuum_toroidal_field.b0

equilibrium_ids.time_slice.global_quantities.psi_axis
equilibrium_ids.time_slice.global_quantities.psi_boundary
equilibrium_ids.time_slice.global_quantities.magnetic_axis.r
equilibrium_ids.time_slice.global_quantities.magnetic_axis.z
equilibrium_ids.time_slice.global_quantities.beta_pol
equilibrium_ids.time_slice.global_quantities.beta_tor
equilibrium_ids.time_slice.global_quantities.beta_tor_norm
equilibrium_ids.time_slice.global_quantities.ip
equilibrium_ids.time_slice.global_quantities.li_3
equilibrium_ids.time_slice.global_quantities.volume
equilibrium_ids.time_slice.global_quantities.area
equilibrium_ids.time_slice.global_quantities.current_centre.r
equilibrium_ids.time_slice.global_quantities.current_centre.z
equilibrium_ids.time_slice.global_quantities.q_axis
equilibrium_ids.time_slice.global_quantities.q_95
equilibrium_ids.time_slice.global_quantities.energy_mhd
equilibrium_ids.time_slice.boundary.type
equilibrium_ids.time_slice.boundary.psi
equilibrium_ids.time_slice.boundary.minor_radius
equilibrium_ids.time_slice.boundary.elongation
equilibrium_ids.time_slice.boundary.triangularity_upper
equilibrium_ids.time_slice.boundary.triangularity_lower
equilibrium_ids.time_slice.boundary.geometric_axis.r
equilibrium_ids.time_slice.boundary.geometric_axis.z
equilibrium_ids.time_slice.boundary.closest_wall_point.r
equilibrium_ids.time_slice.boundary.closest_wall_point.z
equilibrium_ids.time_slice.boundary.closest_wall_point.distance
equilibrium_ids.time_slice.contour_tree.node.critical_type
equilibrium_ids.time_slice.contour_tree.node.r
equilibrium_ids.time_slice.contour_tree.node.z
equilibrium_ids.time_slice.contour_tree.node.psi

equilibrium_ids.time_slice.boundary.outline.r
equilibrium_ids.time_slice.boundary.outline.z

equilibrium_ids.time_slice.profiles_2d.type.index
equilibrium_ids.time_slice.profiles_2d.grid_type.index
equilibrium_ids.time_slice.profiles_2d.grid.dim1
equilibrium_ids.time_slice.profiles_2d.grid.dim2
equilibrium_ids.time_slice.profiles_2d.r
equilibrium_ids.time_slice.profiles_2d.z
equilibrium_ids.time_slice.profiles_2d.psi         #extended to vacuum if using export_field_extension
equilibrium_ids.time_slice.profiles_2d.j_phi
equilibrium_ids.time_slice.profiles_2d.b_field_r   #extended to vacuum if using export_field_extension
equilibrium_ids.time_slice.profiles_2d.b_field_z   #extended to vacuum if using export_field_extension
equilibrium_ids.time_slice.profiles_2d.b_field_phi
```

### `radiation`

Available only for SPI simulations using markers.

```text
The GGD structures are identical to those of plasma_profiles/1, otherwise

radiation_ids.time
radiation_ids.process.ggd.time

radiation_ids.process.identifier.name
radiation_ids.process.identifier%description
radiation_ids.process.identifier%index
radiation.process[0].ggd.ion[0].emissivity.coefficients
```
