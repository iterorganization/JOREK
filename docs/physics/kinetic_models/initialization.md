# Kinetic particle initialization 

With initialization we mean what the particles get as original distribution before the simulation with kinetic_main starts. Initialization can be done in 3 ways, either nothing is initialized, the particles are initialized from a restart file, or the particles are initialized according to some initial distribution.

# No initialization 
Simply set restart_particles=.f. in the namelist, and set nothing for part_group_configs(...)%init_function and %init_pdf (both default to `'none'`). This means that at the start of the simulation with kinetic_main, there will be no simulated particles yet. This is for instance useful for doing neutrals simulations which will be puffed/generated through (wall) recombination automatically during the simulation. 

This only applies to groups with coupling scheme `'ics'` or `'ncs'`, for which initialisation is always skipped. Groups with coupling scheme `'epf'` or `'rep'` **must** set an `init_function` when restart_particles=.f., otherwise kinetic_main aborts (see [below](#initialize-the-particles-according-to-a-given-distribution)).

# Restarting the particles 
To restart the particles from a file, set restart_particles=.t. During initialization, kinetic_main will automatically look for a part_restart.h5 file in the simulation folder and load all the groups from that file directly into the simulation. To do this, you need the corresponding groups in your namelist input file as well. This means you need to have entries in your namelist input with the same %id, %Z, %mass and %type as before the groups before restart, but you can vary each group's %n_particles (although at the moment you can only increase it, not decrease it).

If you want to do your new simulation with more or less part_groups than are present in the restart_file, it is advised to specify your intention by setting the namelist input "part_groups_in_use" to the list of IDs that you want to use in your new simulation (e.g. `part_groups_in_use = "D01" "N01"`). If you want to drop groups from the part_restart.h5 file, you **have** to specify "part_groups_in_use". You will get messages in your output file stating which groups have been loaded and which groups have been dropped.

The part_restart.h5 files are generated automatically at the end of each successful simulation with kinetic_main, and every 500 jorek steps, an intermediate "interim_part_restart.h5" file will be generated (that overwrites a previous interim restart file if that exists, to save on disk space as these restart files can be GBs of data). 

# Initialize the particles according to a given distribution 

When restart_particles=.f., kinetic_main allocates `n_particles` superparticles for each group and then calls `initialise_particles_for_sim` (in `particles/initialisers/mod_particle_group_init.f90`). The initial particle set is written to `part_restart_init.h5`. Set `init_particles_only = .t.` to stop right after that, e.g. to generate a `part_restart.h5` for later runs.

For each group, initialisation is controlled by two namelist settings:

* `part_group_configs(i)%init_function`: selects the sampler, i.e. the velocity-space distribution (and, for `'experimental'`, the real-space distribution as well).
* `part_group_configs(i)%init_pdf`: selects an optional real-space (R,Z) distribution, applied by rejection sampling. It is ignored by `'experimental'`.

What happens depends on the group's coupling scheme:

| coupling_scheme | behaviour |
|-----------------|-----------|
| `ics`, `ncs` | initialisation skipped, particles are born during the simulation (puffing, recombination, ...). `init_function` / `init_pdf` must stay `'none'`. |
| `non` | initialisation skipped with a warning (not fully implemented). |
| `epf`, `rep` | the group is initialised using `init_function` and `init_pdf` as described below. |

## Initialisation functions (`init_function`)

| init_function | particle type | description | required parameters |
|---------------|---------------|-------------|---------------------|
| `'maxwell'` | `particle_kinetic_leapfrog` | Maxwellian energy distribution at constant temperature, uniform pitch. Real-space distribution set by `init_pdf`. Particles are copied `n_phi_planes` times around the torus. | `T_maxwell` [eV], `n_phi_planes` (≥ 1) |
| `'re_gaussian'` | `particle_kinetic_relativistic` | Gaussian energy distribution (monoenergetic if `re_std_energy` is very small) with a single pitch value. Real-space distribution set by `init_pdf`. Charge is set to -1. | `re_energy` [eV], `re_std_energy` [eV], `re_pitch` (= p_par/p_tot, must be non-zero) |
| `'experimental'` | `particle_kinetic_leapfrog` | Samples from a 4D experimental distribution F(R, Z, pitch, E) read from `experimental_dist.h5` in the run directory, see [energetic ions](./energetic_ions/energetic_ions.md#importing-a-fast-particle-experimental-distribution). `init_pdf` is ignored. | `n_phi_planes` (≥ 1) |

`'re_gaussian'` is the only valid choice for `particle_kinetic_relativistic`, and it cannot be used with `particle_kinetic_leapfrog`. Other particle types are currently rejected. All other `init_function` values (including the old `'basic'`, now `'re_gaussian'`) give an error.

For all of the above, `part_group_configs(i)%n_particles_total` (the number of physical particles, i.e. the sum of all weights) **must** be set. All superparticles get the same weight, `n_particles_total / n_particles`. This replaces the former `num_re` parameter for `'rep'` groups.

Missing or inconsistent parameters are caught at start-up by `check_compatibility_and_determine_coupling_schemes` (`particles/mod_coupling_settings.f90`).

## Spatial distributions (`init_pdf`)

The spatial distributions are implemented in `particles/initialisers/mod_spatial_pdf.f90`. Positions are drawn uniformly in the (R,Z) bounding box of the domain and accepted with probability f(R,Z) ∈ [0,1]:

| init_pdf | f(R,Z) | notes |
|----------|--------|-------|
| `'none'` | 1 | uniform in the domain (default) |
| `'RZ'` | $1/(1+r)^2$ | $r$ is the distance to the magnetic axis in metres; broad, peaked on axis |
| `'analytical'` | $\max(1-(r/a)^2, 0)^2$ | $a$ is the minor radius of the LCFS |
| `'itpa_tae'` | fit of the EP profile in $\psi_N$ | reproduces the ITPA TAE benchmark [A. Könies et al 2018 Nucl. Fusion 58 126027] |
| `'current'` | $(j_\phi/R - \min)/(\max - \min)$ | follows the toroidal current density, normalised to its range in the equilibrium. Needs a model with `var_zj` (e.g. model600), not available in full MHD. |

### Example

```
 restart_particles   = .f.
 init_particles_only = .f.

 part_group_configs(1)%id                = 'REs'
 part_group_configs(1)%Z                 = -1
 part_group_configs(1)%coupling_scheme   = 'rep'
 part_group_configs(1)%type              = 'particle_kinetic_relativistic'
 part_group_configs(1)%n_particles       = 1e2
 part_group_configs(1)%mass              = 0.000548579870184805
 part_group_configs(1)%init_function     = 're_gaussian'
 part_group_configs(1)%init_pdf          = 'RZ'
 part_group_configs(1)%n_particles_total = 1.175403e17
 part_group_configs(1)%re_energy         = 9.8e4
 part_group_configs(1)%re_std_energy     = 5e2
 part_group_configs(1)%re_pitch          = -1
```

(taken from the `particle_rep_600` regression test).

## Developer notes

### Code layout

| file | content |
|------|---------|
| `mod_particle_group_init.f90` | `initialise_particles_for_sim` / `initialise_group`: loops over the groups, checks `init_function` against the particle type, builds the spatial pdf, calls the sampler and sets the weights. Replaces `mod_initialise_particles` and `initialisers_RE`. |
| `mod_initialiser_functions.f90` | the samplers: `initialise_particles_maxwell`, `initialise_particles_maxwell_phiplanes`, `initialise_particles_re_gaussian`, `initialise_particles_in_phase_space`, and helpers such as `adjust_particle_weights`. Formerly `initialisers_base`. |
| `mod_spatial_pdf.f90` | the `spatial_pdf` type, the built-in rejection functions and `spatial_pdf_from_name`. |
| `mod_import_experimental_dist.f90` | `initialise_particles_experimental` (formerly `import_particles`). |

Renamed routines: `initialise_particles_H_mu_psi` → `initialise_particles_maxwell`, `import_particles` → `initialise_particles_experimental`.

### The `spatial_pdf` type

The samplers no longer take a rejection function and its variable list as separate arguments (`uniform_space_rej_f`, `uniform_space_rej_vars`). They take one `space_pdf` argument of type `spatial_pdf`, which keeps the two together:

```fortran
type :: spatial_pdf
    procedure(rej_f), nopass, pointer :: f => null()
    integer,          allocatable     :: vars(:)
    logical                           :: needs_grad = .false.
end type spatial_pdf
```

`f` must be a **pure** function with the interface

```fortran
pure function rej_f(n, P, gradP)
    integer,                intent(in) :: n
    real*8, dimension(n),   intent(in) :: P       ! values, ordered as vars(:)
    real*8, dimension(3,n), intent(in) :: gradP   ! d/dR, d/dZ, d/dphi (only if needs_grad)
    real*4                             :: rej_f   ! acceptance probability in [0,1]
end function
```

`vars(:)` lists what is passed in `P`: a positive entry is a JOREK variable index, `-1` = R, `-2` = Z, `-3` = phi (`0` = constant 1). `vars` **must be in ascending order** (geometric entries first), which `check_spatial_pdf` enforces. Gradients are only interpolated when `needs_grad = .t.`.

To use a custom profile in your own program (see `particles/examples/seed_particles.f90`):

```fortran
type(spatial_pdf) :: space_pdf
space_pdf%f    => my_rej_function
space_pdf%vars =  [-2, -1]          ! Z, R
call initialise_particles_maxwell(particles, sim%fields, pcg32_rng(), mass, &
     uniform_space=.true., space_pdf=space_pdf, charge=1)
```

To make a new profile selectable via `init_pdf`, add the rejection function to `mod_spatial_pdf.f90` and add a `case` for it in `spatial_pdf_from_name`. Do any normalisation against the equilibrium there as well (see `calibrate_current_pdf`).
