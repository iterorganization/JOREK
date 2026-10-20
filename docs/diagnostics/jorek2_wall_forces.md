---
title: "jorek2_wall_forces"
nav_order: 32
render_with_liquid: false
parent: "General"
---

# jorek2_wall_forces

This program reads JOREK restart files and calculates the total wall forces. **You need to run your simulations in full free-boundary mode (including the n=0 mode)**. Otherwise you cannot know the vacuum fields outside the plasma domain. For more information on how to run in free-boundary see [here](jorek-starwall-faqs.md). 

It involves the computation of expensive volume integrals in the plasma, so run it in parallel with several MPI processes for large cases. Submit the same job scripts that you would use for JOREK (just replacing "jorek_modelXXX" by "jorek2_wall_forces". To compile it just type

   make -j 8 jorek2_wall_forces

Update 30.11.22 New input file 

To select the given JOREK restart files create a file named wall_forces.nml
with the following format
   &wall_forces
   istart     = 100  ! start index of restart file
   iend       = 1000 ! end index of restart file
   delta_step = 40   ! 
   scale_fact = 1.01 ! scaling factor of wall (1.01 default)
   n_phi_int  = 64   ! number of points for integration in the toroidal direction (to get B-field)
   
   /
The forces are exported to the total_wall_forces_istart..iend..scale_fact.dat file

For older versions the wall_forces.nml just requires one line with istart iend delta_step:
   0 1000 100
The result is stored in total_wall_forces.dat

IMPORTANT: To compute the total wall forces, an integral over a close surface containing the STARWALL wall is performed ([equation 7 in this paper](https://iopscience.iop.org/article/10.1088/1741-4326/aa8876/meta)). The routine assumes a given ordering of the triangles describing this surface, which is an scaled value of the STARWALL wall. Thus depending of how this ordering is done in STARWALL, the triangle normals may flip sign and so the computed forces. TODO: enforce always outward pointing normals in the routine to ensure the correct sign of the total forces
