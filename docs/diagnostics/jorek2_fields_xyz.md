---
title: "jorek2_fields_xyz"
nav_order: 27
render_with_liquid: false
parent: "General"
---

# jorek2_fields_xyz

This program reads JOREK restart files and calculates the magnetic field at arbitrary (x,y,z) points. The particularity with respect to jorek2_postroc is that you can get the fields also outside the JOREK grid. For example, it could be used to produce GEQDSK files from JOREK restart files.

To get the complete fields **you need to run your simulations in full free-boundary mode (including the n=0 mode)**. Otherwise you cannot know the vacuum fields outside the plasma domain. For more information on how to run in free-boundary see [here](jorek-starwall-faqs.md). If case are running with freeboundary=.f. then you only get the fields produced by the plasma currents.

It involves the computation of expensive volume integrals in the plasma, so run it in parallel with several MPI processes for large cases. Submit the same job scripts that you would use for JOREK (just replacing "jorek_modelXXX" by "jorek2_fields_xyz". To compile it just type

   make -j 8 jorek2_fields_xyz

As input you must provide the points in cartesian coordinates in the file xyz.nml in the format

   number_of_points
   x_1       y_1       z_1  
   x_2       y_2       z_2  
    :         :         :   
    :         :         :   
   x_n       y_n       z_n  

DO NOT specify the points exactly at the STARWALL's coil/wall triangles. Otherwise singularities at those points may occur!

To select the given JOREK restart files create a file named steps.nml with 3 integers in the first line (istart, iend, delta_step). For example read the restart files from 100 up to 520 with intervals of 10 steps, just write this first line in steps.nml

   100  520  10

The fields are exported in the file "fields_xyz.dat" and are also separated in coil/wall/plasma contributions. 

Note that to calculate the require plasma fields, n_plane typically should be more than 60 (even for 2D)
