---
title: "Energetic Ion"
nav_order: 7
parent: "Kinetic Particle Module"
layout: default
render_with_liquid: false
---

# Energetic ions

A species of suprathermal ions can be included in JOREK via the kinetic particle framework. To simulate the interaction of energetic particles (EPs) with MHD modes a hybrid kinetic-MHD model is used [1]. The plasma is split into a bulk part, which is treated with an MHD model, and the energetic particles, which are treated kinetically. More details about the JOREK implementation and applications can be found in [2].

## Full-f

JOREK uses a full-f formulation for the EPs. Their distribution function is modeled by $N_\mathrm{p}$ marker particles:

$$
    f(\mathbf{r}, \mathbf{v}, t) = \sum_{i=1}^{N_\mathrm{p}} w_i \, \delta(\mathbf{r}-\mathbf{r}_i(t)) \, \delta(\mathbf{v}-\mathbf{v}_i(t))
$$

The markers are evolved in time with one of the available pushers, e.g. the Boris pusher.

## Coupling scheme

The coupling between EPs and bulk MHD fluid is provided by a pressure coupling scheme. The EP pressure tensor $\boldsymbol{P}_\mathrm{h}$ enters the momentum balance equation according to

$$
  \rho \frac{\partial \mathbf{u}}{\partial t}=-\rho \mathbf{u} \cdot \nabla \mathbf{u}+\mathbf{J} \times \mathbf{B}-\nabla p +  \left[ \nabla\cdot\underline{\boldsymbol{P}_\mathrm{h}}\right]_\perp
$$

Instead of the full tensor, usually the CGL form

$$
  \underline{\boldsymbol{P}}_\mathrm{h} =  \left[  P_{\parallel} \mathbf{b}\mathbf{b} + P_{\perp}\left( \underline{1}-\mathbf{b}\mathbf{b} \right)  \right]
$$

is used with components

$$
  P_\parallel = m \int \mathrm{d}^3 \mathbf{v} \ f(\mathbf{r}, \mathbf{v}, t) v_\parallel^2
$$

$$
  P_\perp = \frac{m}{2} \int \mathrm{d}^3 \mathbf{v} \ f(\mathbf{r}, \mathbf{v}, t) v_\perp^2
$$

In reduced MHD, it is also possible to switch to a scalar pressure for the EPs by setting `use_pcs_full = .false.` and `use_pcs = .true.` In the full MHD models, only the CGL form is implemented and always switched on by `use_pcs = .true.`

The pressure tensor is provided by the six components

$$
  \boldsymbol{P}_\mathrm{h}^{RR}, \ \boldsymbol{P}_\mathrm{h}^{ZZ}, \ \boldsymbol{P}_\mathrm{h}^{\phi\phi} \ \boldsymbol{P}_\mathrm{h}^{\phi\phi}, \ \boldsymbol{P}_\mathrm{h}^{RZ}, \ \boldsymbol{P}_\mathrm{h}^{R\phi}, \ \boldsymbol{P}_\mathrm{h}^{Z\phi}
$$

which are then used in the MHD system of equations in `mod_elt_matrix_fft`.

## Example

An example use case can be found in `particles/examples/tae_phase_space_project.f90` with the input in `namelist/model307/itpa(_eq)` or `namelist/model710/itpa_full(_eq)`. This is the setup for the ITPA EP benchmark case producing a TAE [3].

Running a simulation including EPs in JOREK works in two steps:

1. Creating the MHD equilibrium using a desired jorek model (`jorek_modelXXX`).  Run the simulation for 0 timesteps. This will produce the file `jorek_restart.h5`.
2. Running the coupled particle simulation. This can be either done with a dedicated particle exectuable (such as `tae_phase_space_project`) or using the generalized program `kinetic_main`.



## Diagnostics

Phase space projections provide very useful and flexible diagnostics for EPs. In the file `particles/diagnostics/mod_phase_space_project.f90` the interface and the functions for phase space projections are extensively commented.

The basic idea of a projection is to calculate any particle quantity as a function of any real or velocity space coordinate. An example is the energy exchange during some specific time period as a function of the parallel velocity and magnetic moment:

$$
  G(v_\parallel, \mu) = \sum_{i=1}^{N_\mathrm{p}} \Delta E_i \ w_i \delta(v_\parallel - v_{\parallel,i}) \delta(\mu - \mu_i)
$$

where $\Delta E_i$ is the energy change of the $i$th marker particle.

In JOREK, these projections can be done in several ways. It is possible to use the function `new_phase_space_projection` and provide the projection function beforehand, or, fill up a value array with projection values for each single marker particle using the subroutine `project_single_particle_x`.




## Initialization of energetic particles

### Consistency with the MHD equilibrium

EPs can contribute to a large fraction of the total pressure. For this reason, in nonlinear simulations it is important to approximately satisfy the condition

$$
  p_\mathrm{eq} = p_\mathrm{EP} + p_\mathrm{bulk}
$$

where $p_\mathrm{eq}$ is the pressure that has been used to compute the MHD equilibrium. This can be achieved for example by using different profiles for the temperature and / or density in the calculation of the equilibrium and in the coupled particle simulation.

### Desired input profiles

EPs can have large orbit widths and do not necessarily remain on a specific flux surface. Therefore, it is not trivial to initialize a distribution of marker particles based on desired profiles as a function of the normalized flux in the full-f scheme. A naive initialization will result in a relaxation of the profiles.

### Importing a fast particle experimental distribution

The option to import an experimental distribution function for fast particles is located in the file `particles/initialisers/mod_import_experimental_dist.f90` (subroutine `initialise_particles_experimental`, formerly `import_particles`).
In kinetic_main it is selected per group with `part_group_configs(i)%init_function = 'experimental'` (with `restart_particles = .f.`). This needs `type = 'particle_kinetic_leapfrog'`, `n_phi_planes` ≥ 1 and `n_particles_total`, and the file must be called `experimental_dist.h5` and be in the run directory. `init_pdf` is ignored, since the spatial distribution is part of the file. From kinetic_main all particles are copied over the $\phi$ planes (`fraction_phi_planes = 1`). See also [Kinetic particle initialization](../initialization.md).
For now, the only option is to use a 4D distribution function, as it is assumed the other two dimensions (gyrophase and toroidal angle) are uniformly distributed. The dimensions then are $(R,Z, \lambda,E)$ where $\lambda$ is the pitch angle defined as  $\lambda = \frac{v_\parallel}{|\mathbf{v}|}$. This distribution function is then 4D rejection sampled, with linear interpolation on the 4D grids. The gyrophase and toroidal angle are then sampled uniformly, but with the option to copy the particles $n_\phi$ amount of times with $\phi$ increasing by $\frac{2\pi}{n_\phi}$ (eliminating toroidal harmonics). This can also be done for a subset of the particles (as to keep some noise perturbation large enough).

The expected HDF5 is structured as follows:
   * 4 1D grids: 'E_1D', 'Pitch_1D', 'R_1D' and 'Z_1D'
   * 1 4D normalized distribution function 'F0_norm'. This cannot have values > 1. Values < 0 are ignored (which doesn't matter for rejection sampling).

The units are:
   * R [m]
   * Z [m]
   * Pitch {-1.0, ...,1.0}. The sign can be somewhat confusing, as sometimes it is given with respect to J and sometimes to B. In the current implementation, the usecase is for AUG, which gives it with respect to J while J is opposite to B, such that $v_\parallel = - \lambda |\mathbf{v}|$ ($v_\parallel$ is in JOREK with respect to B). It is trivial to change this in the code (or in the HDF5). It is advised to check with diagnostics whether you took the correct sign.
   * E [eV]


## References:
[1] W. Park et al., Phys. Fluids B 4, 2033 (1992).

[2] T. Bogaarts et al., Phys. plasmas 29, 122501 (2022).

[3] A. Könies et al., Nucl. Fusion 58, 126027 (2018).






