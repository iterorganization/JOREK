# Documentation JOREK

## Hot tail tracker

This section explains the use of the "Hot tail tracker" framework, which aims to
simulate the primary Runaway Electron (RE) generation during the Thermal Quench
(TQ) phase of a disruption through the tracking of test electron markers. This
tool, presented in [1], primarily relies on:

- **Guiding Centers (GCs) relativistic model**: the gyrophase is neglected
  `>> class(particle_gc_relativistic)`
- **Test particle approximation**: markers are influenced by MHD fields, but not
  the opposite
- **RK4 pusher [3][4]**: pushes markers in time according to electromagnetic
  fields, with the radiation reaction force included
  `>> runge_kutta_fixed_dt_gc_push_jorek_radreact`
- **Collision operator [5]**: accounts for small-angle Coulomb collisions
  `>> ccoll_gc_relativistic_push`
- **Reweighting method [1]**: Correlated Russian Roulette method to maintain a
  satisfactory statistic at high energy
  `>> correlated_russian_roulette`
- **Recycling method [2]**: recycles lost markers
  `>> recycle_lost`

---

## Compilation and run

For now, only the GNU openmp compiling is possible (the rejection-sampling
initialization does not work with the INTEL oneapi one). See the Compiling and
Running section for more information. The compilation is made by:

```
>> make -j 8 hot_tail_tracker
```

**Known issue**: error during compilation `Error: no such instruction: 'vmovw ...'`
is solved by deleting `-march=native` in `Makefile.inc`.

The launch must be performed in the folder in which there are:

- `jorek*.h5` files for the fields import
- `ccolldata`, the pre-calculated look-up table for collisions, obtained by
  compiling and running `particles/util/ccoll_generate_L0L1.f90`
- the input file (for instance `jorek_namelist`)

The simulation is performed by `particles/examples/hot_tail_tracker.f90`. It can
be run using:

```
>> mpirun ./hot_tail_tracker < jorek_namelist > jorek_output_hot_tail_tracker
```

The simulation is hybrid OpenMP/MPI parallelized.

**WARNING**: the number of markers simulated must be divisible by the number of
MPI processes.

**Known issue**: as for other simulations with kinetic particles, the use of many
particles may lead to memory issues. The stack size should be increased by adding
the following to the jobscript:

```
>> ulimit -s unlimited
>> export OMP_STACKSIZE=1024m
```

**Thread binding**: under Open MPI a rank pinned to a single core shares that
core between all its OpenMP threads, which silently removes the OpenMP speedup.
Check the binding of the jobscript if the threads do not scale.

---

## Output

During the simulation, 5 different individual marker data files are created:

- `jorek_output_hot_tail_tracker`: the output logfile
- `diag.h5`: the individual state of markers at different times, with notably
  `p` in units of `m_e c`, the norm of the momentum, `xi` the cosine of the pitch
  angle `xi = cos(theta) = p_par/p`, and `(R, Z, Phi)` the JOREK coordinates
- `Jcoll_....h5`: the field seen by each marker (density, temperature, electric
  field) at different times. Stored in single precision to halve the file size
- `Jlosses_....h5`: the list of the last state of each lost marker
- `part_restart.h5`: the file for restarting the kinetic simulation

The suffix of `Jcoll_...` and `Jlosses_...` is `K{K1}_{K2}_np{nprt}_dt{timestep}ns`,
with `[K1, K2]` the energy window in which markers are initialized, `nprt` the
number of markers used, and `timestep` the particle time step.

Two points when reading the loss file:

- `iprt` is the **global** marker index, so loss records are comparable across
  MPI processes.
- Only physical losses are recorded. A marker discarded by the momentum jump
  detector is recycled without being written to the loss file, so that the loss
  diagnostic is not polluted by numerical artefacts. Those are counted separately
  in the logfile.

The logfile ends with a per-process summary:

```
Wall time                                   <- the figure for the run duration
CPU time (all threads)                      <- summed over OpenMP threads, so larger
Nb jumps / Weight jumps                     <- pushes rejected by the jump detector
Nb invalid markers found at push entry      <- markers already outside the mesh
Nb failed recycling attempts
Nb stale marker indices after field reloads
```

Three warnings may appear during the run:

| Message | Meaning |
|---|---|
| `failed recycling attempts so far` | No donor was found for a dead marker. A few are harmless. A growing count means the donor conditions are too strict: relax `p_cutoff_recycle` or enable `use_p_far_sweep` |
| `the grid changed at the field reload` | The JOREK mesh differs between two restart files. The affected markers were invalidated and recycled rather than reading out of bounds. Expected if the MHD run was refined |
| `loss buffer reached ... flushing early` | Normal during a quench. Raise `max_loss_per_rank` if memory allows and fewer writes are wanted |

---

## Parameters

The parameters of the simulation can be directly changed in
`hot_tail_tracker.f90`.

### Simulation set up

- **`dt_marker`**: the time step to push particles, in seconds. 1 to 10 ns is a
  good start
- **`duration`**: duration of the simulation, in seconds
- **`nprt`**: number of markers to track. 1e4 to 1e5 is a good start. It must be
  divisible by the number of MPI processes
- **`time_startsim`**: time to start the simulation, not accounted for if
  `restart_particles = .TRUE.`
- **`nstep_write`**: number of times that the diagnostics are written
- **`nstep_reweight`**: number of times that a reweighting step is performed
- **`average_coll_fields`**: flag to set whether the output density, temperature
  and electric field seen by markers is averaged between two written steps
  (`.TRUE.`) or not (`.FALSE.`)
- **`restart_particles`**: flag to indicate whether the simulation starts from a
  `part_restart.h5` file (`.TRUE.`) or not (`.FALSE.`)
- **`iseed`**: `0` seeds the random number generator from the OS entropy, so
  results differ from run to run. A positive value gives a reproducible per-rank
  seed. Note that gfortran derives the per-thread states from the master state,
  so a reproducible run also requires a fixed `OMP_NUM_THREADS`

### Initialization and reweighting parameters

- **`Rbound`**: window of the R coordinate within which markers are initialized.
  Clipped to the mesh bounding box at runtime
- **`Zbound`**: window of the Z coordinate within which markers are initialized.
  Clipped likewise
- **`Phibound`**: window of the Phi coordinate within which markers are
  initialized
- **`Ekinbound`**: window of the kinetic energy, in eV, within which markers are
  initialized
- **`Pitchbound`**: window of the pitch angle, in radians, within which markers
  are initialized. `0` is aligned with the magnetic field and `pi` is
  anti-aligned. The angle is sampled uniformly in its cosine
- **`total_weight`**: total weight of markers. Advice: set it to the total number
  of electrons in your plasma if you simulate all the phase space
- **`range_weight`**: range of weight to scan. Markers with a weight below
  `total_weight/range_weight` are not initialized. The same threshold also
  defines a dead marker during the run, and the minimum weight a marker must have
  to be used as a recycling donor. It therefore bounds the geometric decay of the
  weights under repeated splitting
- **`sigma_psi`**: standard deviation of the PI function,
  `PI(psi) = exp(-(psi-psi_axis)^2/(2*sigma_psi^2))`, allowing to place more
  markers spatially where PI is high. Set a high value, for instance 1000, to
  have a uniform initialization
- **`jump_coll`**: the collision operator diverges at low momentum. It may happen
  that the momentum of a particle jumps to an unrealistically high value after a
  collision. This parameter flags a jump when `p(t)` is above the operator
  threshold `p_th` and `p(t+dt) > jump_coll * p_th`, and discards the marker
  otherwise. 3 is a reasonable value
- **`E_par_max`**: upper bound assumed on `|E_par|`, in V/m, used to build the
  jump detection threshold of the RK4 push itself,
  `delta_jump = e * E_par_max * dt_marker / (m c)`. A push changing `p/(mc)` by
  more than that is treated as a numerical artefact. Raise it if the case
  genuinely carries a stronger parallel field, otherwise real acceleration will
  be discarded as a jump
- **`p_cutoff`**: momentum boundary for the reweighting method (`p_b` in [1]). If
  a marker crosses this limit, it plays the Russian roulette. Set a negative
  value to not use the reweighting method
- **`p_far`**: momentum limit for the reweighting method (`p_f` in [1]) over which
  duplication is allowed after a kill by the Russian roulette
- **`k`**: k for the reweighting method
- **`l`**: l for the reweighting method, `proba_save = l/k`
- **`age_split`**: number of reweighting processes that a marker should wait after
  its birth to be a candidate for duplication
- **`use_psi_res`**: `.TRUE.` activates the condition of duplication close to
  `psi_res` from the death of the killed marker. May help to avoid clustering
- **`use_p_far_sweep`**: `.TRUE.` allows `p_far` to be reduced by `delta_p_far` if
  no candidate is found after `nprt` draws. Activate it if finding a duplication
  candidate is too hard or impossible
- **`psi_res`**: if `use_psi_res = .TRUE.`, defines the psi distance for the
  duplication condition
- **`delta_psi_res`**: if no duplication candidate is found after `nprt` draws,
  `psi_res = psi_res + delta_psi_res`
- **`delta_p_far`**: if `use_p_far_sweep = .TRUE.` and no duplication candidate is
  found after `nprt` draws, `p_far = p_far - delta_p_far`
- **`p_cutoff_recycle`**: cutoff momentum for the recycling of losses, over which
  a marker can be duplicated

### Memory

- **`max_loss_per_rank`**: maximum number of buffered loss records per process
  before an early flush is forced. `write_loss_to_hdf5` gathers every process
  buffer on rank 0, so the root allocation is
  `n_mpi * max_loss_per_rank * 8 arrays * 8 bytes`, that is 0.41 GB at 32
  processes with the default of 200000. Without a cap the buffers grow with the
  loss rate times the number of steps between two scheduled writes, which reaches
  several GB during a thermal quench

---

## Analysis

A set of python codes is available for the output data analysis: git@github.com:LouisP-runaway/analysis_hot_tail_tracker.git
The user chooses the analysis via the `main.py` script.

---

## References

[1] Puel, *Kinetic modelling of hot tail runaway electron generation during
plasma disruptions using the JOREK code*, NF (2025)

[2] Puel, NF (2026)

[3] Cristian Sommariva, NF (2017)

[4] Xin Tao, PoP (2007)

[5] Konsta Särkimäki, *Computer Physics Communications* (2018)

The header of `hot_tail_tracker.f90` cites two DOIs,
[10.1088/1741-4326/ae2d70](https://doi.org/10.1088/1741-4326/ae2d70) for the
implementation and verification, and
[10.1088/1741-4326/aea175](https://doi.org/10.1088/1741-4326/aea175) for the
ITER baseline H-mode 15 MA case the shipped parameters target.
