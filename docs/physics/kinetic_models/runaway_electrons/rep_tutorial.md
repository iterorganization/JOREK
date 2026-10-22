---
title: "Running with Kinetic Runaway Electrons"
nav_order: 8
parent: "Kinetic Runaway Electrons"
layout: default
render_with_liquid: false
---

# Getting Started with full-kinetic Runaway Electrons on kinetic_main (WIP 2026 Feb, Chizhou) 

This work is based on kinetic_develop branch [Introduction to the kinetic framework](./../particle_introduction.md), where RE super particles are evolved in full orbits and coupled with the fluid plasma. The compiling of the binary kinetic_main can be found in [Introduction to the kinetic framework](./../particle_introduction.md), which is similar to jorek2_main. The Strumpack solver have been usable for kinetic simulations. More general introduction to kinetic models in JOREK is given in [particles](./../particle.md)

# Initialize the REs

Apart from the jorek_restart.h5, an part_restart.h5 file is required by any simulations with particles. Without particle restart, kinetic_main runs the same way as jorek2_main. The guide for initialization can be found in [initialise_particles_in_phase_space](./../initialization.md). To start with, ''test_initialisation_phase_space'' in particles/examples is a good example for initializing REs, where the energy, pitch, and toroidal coordinate spectrums are initialized uniformly within a specified range, while the poloidal distribution can be set to follow the current density in a fluid restart file ''jorek_equilibrium_restart.h5''.

**ATTENTION:** the example scripts are now outdated. The ''group%id'' should be three letters in latest kinetic_main, but it's "001" in the script. You will need to edit the group id in the particle restart file. 

# Run main_kinetic

''main_kinetic'' works similarly as ''jorek2_main'', but it can **ONLY** restart from a ''jorek_restart.h5'' **AND** a ''part_restart.h5'' instead of creating a new equilibrium with ''restart=.f.''. Apart from the particle restart file, You need to add some extra information about the particle groups in namelist (also see [Introduction to the kinetic framework](./../particle_introduction.md)). An RE-specified example is given below. You can add it to the end of your fluid namelist:

```
 restart_particles = .true.                                 ! If false, run as pure fluid without REs
 ! Particle properties in new format
 part_group_configs(1)%id                          = 'REs'  ! group id you give when initializing particles, MUST be 3 letters
 part_group_configs(1)%Z                           = -1
 part_group_configs(1)%coupling_scheme             = 'rep'  
 part_group_configs(1)%type                        = 'particle_kinetic_relativistic'
 part_group_configs(1)%n_particles                 = 1e6    ! number of super particles
 part_group_configs(1)%mass                        = 0.000548579870184805     
 part_group_configs(1)%n_particles_total           = 1.175403e16      ! number of physical REs, doesn't work if initialized according to current
 part_group_configs(1)%use_kin_recombination       = .false.

! Projecting (particle field -> MHD) smoothing
 filter_hyper        = 1.d-9
 filter_hyper_n0  = 1.d-9
 filter_par            = 5.d-1
 filter_par_n0       = 5.d-1
```

When bringing in particles, part of ohmic current is automatically converted to RE current carried by particles by keeping a consistent total current [runaway_electrons](./runaway_electrons.md). This might lead to noise so it is recommended to run the simulation for some time without bringing in new physics to smooth the fluid model. To create the scenario where current is carried purely by REs, you may give a high value to ''eta'' and set ''E=0'' on line 199 of ''particles/mod_fields.f90'', so that the ohmic current decays while the RE particles won't be accelerated by the artificially induced electric field.   



   



