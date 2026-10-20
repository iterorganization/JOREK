---
title: "jorek2_fast_camera"
nav_order: 25
render_with_liquid: false
parent: "General"
---

## jorek2_fast_camera

### Parallelisation
This diagnostic is parallelised in the number of lines-of-sight (ie. the number of pixels).  
Hence to run it, you need to use a lot of MPI processes (not threads). 

### ADAS data
There is one file needed to run the diagnostic: my_pec.dat, which needs to be in the running directory.  
This file contains the ADAS emission for the Dalpha lines, which is determined by the electron density and temperature.  
This is one of the main improvement that can be brought to this code: to use the exact ADAS data, for Dalpha, but potentially for impurity radiation as well.  
Otherwise, just the jorek_restart.h5 file is needed.

### Camera view and geometry
The diagnostic is currently set up for the MAST fast camera view.  
This can be changed inside the file itself (`jorek/diagnostics/jorek2_fast_camera.f90`), with the following parameters:  
`pixel_dim` : size of pixels  
`n_pixels_hor` : number of horizontal pixels  
`n_pixels_ver` : number of vertical pixels  
`X_cam` : camera position  
`Y_cam` : camera position  
`Z_cam` : camera position  
`focus` : focal length of lens  
`toroidal_angle_location` : toroidal angle position (easier than using X,Y,Z, so that Z can always be zero)  
`angle_hor` : orientation of camera (horizontal)  
`angle_ver` : orientation of camera (vertical)  
`light_saturation` : factor to change maximum saturation of camera (otherwise highest pixel light will be used)  
`colormap` : normally just Black-and-white, but can also use heatmap or rainbow...  
''''  
  

All input parameters can be set either in the code `jorek2_fast_camera.f90` before compilation, or at run-time, in the input file `fast_cam.nml`, which must be placed in the running directory.  
Examples of these files for JET and MAST can be found in `jorek/util/fast_camera/`



### Neutrals
The Dalpha emission uses the neutrals density. If model500 is used, it uses the actual rho_n variable for this, but if we are using model303 for example, the neutrals density is set by hand, inside the code (or inside the `fast_cam.nml` input file):  
`  rhon_prof(1) = 1.d-3  ! Core value [10^20 m^(-3)]`  
`  rhon_prof(2) = 1.d-1  ! SOL value [10^20 m^(-3)]`  
`  rhon_prof(3) = 3.d-1  ! Lower-divertor value [10^20 m^(-3)]`  
`  rhon_prof(4) = 3.d-1  ! Upper-divertor value [10^20 m^(-3)]`  
`  rhon_prof(5) = 0.10   ! Edge-psi tanh width [psi_norm]`  
`  rhon_prof(6) = 0.95   ! Edge-psi tanh position [psi_norm]`  
`  rhon_prof(7) = 0.20   ! Lower-divertor tanh width [m]`  
`  rhon_prof(8) = -0.05  ! Lower-divertor tanh position relative to Z_xpoint(1) [m]`  
`  rhon_prof(9) = 0.20   ! Upper-divertor tanh width [m]`  
`  rhon_prof(10)= +0.05  ! Upper-divertor tanh position relative to Z_xpoint(2) [m]`  


### Line-of-sight integration
The integration along a line of sight is done with the step specified by the user (in meters), eg.  
`step   = 1.d-3`  
But it is possible to use a varying integration step-length, which depends on whether we are entering the plasma core (low emission on MAST), or whether we are at the plasma edge (at the moment `psi_n=0.6`).  
This is completely arbitrary, and can be switched on by the user, by changing the variable  
`go_faster_in_core = .true.`  


### Output format
At the moment, the output is written simply as a .ppm image, which is a very simple format to write.  
It may be desirable to use a .vtk output format as well, or maybe h5...  


### Reflections
At present reflections are included only for a toroidally symmetric wall.  
This can be done by specifying the limiter contour in the JOREK input file:  
`n_limiter`  
`R_limiter(:)`  
`Z_limiter(:)`  
and by setting the reflections flag and reflective coefficient in the code (or the `fast_cam.nml` input file):  
`include_reflections`  
`reflective_coefficient`  
Only the first reflection is taken into account.  
Developing the code to include reflection on CAD structures would require more work, but the line-of-sight structuration in the code is now ready to do this relatively easily.


### Example
Here are examples of the fast camera on MAST, and the wide-angle view on JET, including reflections  
![](fast_camera_mast2.png?400)  
![](fast_camera_jet2.png?400)
