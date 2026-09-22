---
title: "Run Stellarator Simulations"
nav_order: 1
parent: "Stellarator"
grand_parent: "Howto"
layout: default
render_with_liquid: false
---

# Setting up a basic stellarator simulation

This wiki page outlines how to run a W7-A stellarator simulation using JOREK, as an introduction to the stellarator models.

The process of setting up a stellarator simulation involves several steps, but the main idea can be summarized as follows:

  - Calculate an equilibrium in VMEC
  - Refine the equilibrium in GVEC
  - Extract the necessary data for JOREK:
    - Run the GVEC to JOREK convertor to get the necessary equilibrium data for JOREK.
    - (If Dommaschk potentials are used, compute this representation of the vacuum magnetic field)
  - Reconstruct the equilibrium in JOREK using model 180
  - Run the simulation in time using model 183

The following sections describe each step in detail using the Wendelstein 7-A stellarator as an example.

## 1. Calculating a stellarator equilibrium in VMEC
Unlike the tokamak models, JOREK does not have a built-in equilibrium solver for stellarators. For this reason, JOREK uses an ideal MHD equilibrium solution with nested as its initial condition. This is computed using VMEC and GVEC. In the future, it should become possible to skip the first step and go directly to (2), generating the equilibrium directly in GVEC.

VMEC2000 is compatible with GVEC. A description of all parameters in VMEC2000/STELLOPT VMEC is available [here](https://princetonuniversity.github.io/STELLOPT/VMEC%20Input%20Namelist%20(v8.47)), alongside more detailed tutorials on how to run the code.

In order to run VMEC, one needs to specify some input parameters in a namelist file, which is then passed to VMEC. A sample namelist file for a low beta (1 Pa on axis) Wendelstein 7-A equilibrium is provided below:

```fortran
&INDATA
lfreeb        = .false.
MGRID_FILE    = 'none'
NFP           = 5,
mpol          = 12,
ntor          = 10,
ntheta        = 60,
nzeta         = 50,
NS_ARRAY      = 12, 47, 133, 201, 501, 1001, 3001 
FTOL_ARRAY    = 1.e-06, 1.e-09, 5.e-10, 1.e-12, 1e-12, 1e-12, 1e-12
DELT          = 0.2,
NITER         = 100000,
NSTEP         = 100,
AM            = 1.0, -0.99
AC            = 3.0 -6.0 3.0 0.0 
NCURR         = 1,
CURTOR        = -17500.0
GAMMA         = 0.0,
PHIEDGE       = 0.08,
raxis      = 1.99,
zaxis      = 0.,
rbc(0,0)=1.99090000E+00,zbs(0,0)=0.00000000E+00,
rbc(1,0)=5.17790000E-04,zbs(1,0)=-5.12860000E-04,
rbc(-2,1)=-1.20200000E-04,zbs(-2,1)=-1.15270000E-04,
rbc(-1,1)=-5.25220000E-04,zbs(-1,1)=-5.75380000E-04,
rbc(0,1)=1.05000000E-01,zbs(0,1)=1.04000000E-01,
rbc(1,1)=-1.3239E-02,    zbs(1,1)=1.388E-02,
rbc(2,1)=1.13000000E-05,zbs(2,1)=2.31320000E-04,
rbc(1,3)=5.61690000E-05,zbs(1,3)=1.04540000E-04,
rbc(2,3)=3.35470000E-04,zbs(2,3)=-3.36410000E-04,
rbc(1,4)=-1.38590000E-04,zbs(1,4)=-1.35450000E-04,
rbc(3,5)=-1.74710000E-04,zbs(3,5)=1.71180000E-04,
/
&end
```

Pressure, along with one of either toroidal current or rotational transform must be specified. The parameter `NCURR` signals which profile will be provided: if `NCURR = 1` then the toroidal current is specified, and if `NCURR = 0` then the rotational transform is specified. By default, the profiles are given as polynomials:
$$f(s) = \sum_n a_n s^n,$$
where $f$ is the profile in question, $s$ is the normalized toroidal flux and the coefficients $a_n$ are provided by the array `AM` for pressure, `AC` for toroidal current and `AI` for iota (not used here). The total toroidal current is given by `CURTOR`.

The parameters `raxis` and `zaxis` specify the cosine modes of the $R$ coordinate of the magnetic axis, and the sine modes of the $z$ coordinate of the magnetic axis (both as functions of toroidal angle), respectively.

The most important parameters used above that were not mentioned yet can be summarized as follows. The parameter `lfreeb` is set to false meaning that a fixed-boundary equilibrium will be calculated, and so no `MGRID_FILE` is necessary. The parameter `NFP` sets the number of periods, while `mpol` and `ntor` specify the number of poloidal and toroidal Fourier modes, respectively, within one period, and `ntheta` and `nzeta` set the poloidal and toroidal resolution, respectively. The radial resolution is given by the array `NS_ARRAY`. VMEC starts by calculating a low-radial-resolution equilibrium, adding more flux surfaces, minimizing energy to get a higher-resolution equilibrium, then adding more flux surfaces, etc. The number of elements in `NS_ARRAY` specifies the number of times the radial resolution is refined, with each element giving the total number of flux surfaces in the corresponding refinement step. In addition, each refinement can be calculated to a different tolerance; each element of `FTOL_ARRAY` gives the tolerance of the corresponding refinement. Finally, for a given toroidal angle $\zeta$ and poloidal angle $\theta$, the $(R,z)$ coordinates of the boundary are specified as a series of modes:
$$R = \sum_{m,n} a_{mn}\cos(m\theta - n\zeta),\quad z = \sum_{m,n} b_{mn}\sin(m\theta - n\zeta).$$
The arrays `rbc(n,m)` and `zbs(n,m)` specify the coefficients $a_{mn}$ and $b_{mn}$, respectively.

Once the parameters have been set in the namelist file, VMEC can be run as follows:
```bash
./xvmec < input.l2stel_w7a_like
```
where `input.l2stel_w7a_like` is the name of the namelist file above.
## 2. Refining the equilibrium in GVEC

### General information
At the current state, GVEC is a fixed-boundary 3D equilibrium code, using the same strategy as VMEC to find the equilibrium. Stand-alone, it can be run with a given boundary shape and a iota and pressure profile given as polynomial expansion. The option to specify the toroidal current profile instead of iota, is not in a workable state yet.

Alternatively, it can use a given VMEC file to extract that information (iota, pressure, total toroidal flux and fourier modes of the boundary). This was done for the Wendelstein 7-A equilibria, since they were computed with a given toroidal current profile.

The unknowns that are found in VMEC are the geometry of the flux surfaces with label $s$, $R(s,\theta,\zeta),Z(s,\theta,\zeta)$ and the potential $\lambda(s,\theta,\zeta)$, which results in the PEST straight-field line poloidal angle $\theta^*=\theta+\lambda$. 

In GVEC, the unknowns are $(X_1,X_2,\lambda)$. This is because the toroidal map $(R,Z,\phi)->(x,y,z)$ was generalized to arbitrary other analytical maps. In the default case of the toroidal map, the unknowns are the same as in VMEC, so $R=X_1,Z=X_2$ and $\lambda$.

Together with the toroidal $\Phi$ and poloidal fluxes $\chi$, the magnetic field is then defined as
$$ B = \Phi'(s)\nabla s \times \nabla(\theta+\lambda) - \chi'(s)\nabla s \times \nabla \zeta 
     = \Phi'(s)(1+\partial_\theta(\lambda)) \nabla s \times \nabla \theta + (\chi'(s)-\Phi'(s)\partial_\zeta(\lambda)) \nabla \zeta \times \nabla s$$
where typically, the toroidal flux is used as a radial variable (in VMEC its the normalized toroidal flux, in GVEC its the square root of the normalized toroidal flux). Then the polodial flux derivative is given by the iota profile, $\chi'(s) = \iota(s)\Phi'(s)$.

### Input parameters

A sample parameter file for the Wendelstein 7-A case is shown here:
```ini
!===================================================!
! initial solution from existing equilibrium
!===================================================!
  whichInitEquilibrium = 1 ! 0: only from input parameters, 1: from VMEC file

  ProjectName = GVEC_W7A_UNSTABLE_1PA_mpol12_ntor10_r41_vmecBC_accGD
    vmecwoutFile=wout_l2stel_w7a_like.nc
    VMECwoutfile_format = 0 !0: netcdf, 1: nemec

  init_fromBConly =F   ! true: only use axis and boundary for X1,X2 (True is default) 
                       ! false: only needed if VMEC data is used, interpolate on full mesh s=0...
  reinit_BC       = -1 ! =-1, keep vmec axis and boundary,
                       ! = 0: keep vmec boundary, overwrite axis (from parameterfile)
                       ! = 1: keep vmec axis, overwrite boundary
                       ! = 2: overwrite axis and boundary 
  init_LA         = F  ! = T: compute lambda from initial mapping (default)
                       ! = F: lambda=0 / from VMEC at initialization 
!===================================================!
! discretization parameters
!===================================================!
  sgrid_nElems=41          ! grid in radial direction (for B-splines)
  sgrid_grid_type=4        ! spacing: 0: uniform, 1: sqrt (finer at edge), 2: s^2 (finer at axis), 
                           ! 3: bump (fine at edge and center), 4: uniform at axis, finer at edge
  degGP  = 7               ! number of gauss points per radial element
  fac_nyq = 4              ! factor to compute the number of fourier interpolation points: 
                           ! (m_nyq,n_nyq) = fac_nyq * (m_max,n_max)
  X1X2_deg  = 5            ! polynomial degree in radial discretization for X1 and X2 variable
  LA_deg    = 4            ! polynomial degree in radial discretization for Lambda variable
  LA_continuity = 3        ! set to LA_deg-1, so that a spline in radial direction is used
  
!!!these parameters can be set, but they are commented out, so they are taken from the vmec file! 
!!!  X1_mn_max = (/ 12, 10/)   ! maximum mumber of fourier modes (mmax,nmax)
!!!  X2_mn_max = (/ 12, 10/)   ! 
!!!  X1_sin_cos = _cos_        ! which fourier modes: can be either _SIN_,_COS_, or _sin_cos_, 
!!!  X2_sin_cos = _sin_        !
!!!  LA_mn_max = (/ 12, 10/)   ! maximum mumber of fourier modes (mmax,nmax)
!!!  LA_sin_cos = _sin_        !
!===================================================!
! Minimization parameters
!===================================================!
  MaxIter =100000             ! maximum number of iterations
  outputIter = 10000          ! number of iterations after which a state & visu is written
  logIter = 1000              ! write log of iterations to screen 
  start_dt=0.5                ! step length of gradient descent (can be reduced by the algorithm)
  minimize_tol=1.0e-12        ! absolute tolerance |Force|<abstol
  PrecondType = 1             ! -1: off(default), 1: on
  MinimizerType = 10          ! 0: gradient descent(default), 10:  new accelerated gradient descent
```
The parameter `vmecwoutFile` gives the path to the VMEC output file from the previous step, which must be read by GVEC. The parameter `ProjectName` gives a prefix for the GVEC output file names, to which `_State_0000_<stepnumber>.dat` is appended to make the output file name. In principle, only the last state file is needed further on.

The gradient descent/accelerated gradient descent (`MinimizerType=0/10`) minimizes the total MHD energy under the given constraints. The gradient descent is still the more robust (but slower) option. The preconditioner (`PrecondType=1`) greatly helps convergence, and allows a step length of the order of 1 (`start_dt`). The abort criterion using the tolerance `minimize_tol` takes the norm of the gradient vector used in the gradient descent. But this does not exactly correspond to the strong force balance, since it depends on the resolution. Since it is not clear what value is 'small enough', here it is set to a very small value, such that GVEC never achieves it and simply runs for `MaxIter` iterations.

For the discretization, one can specify the number of elements in radial direction (`sgrid_nElems`), which is then used for the B-Spline of certain degree (`X1X2_deg` for the $X_1,X_2$ variables,`LA_deg` for the $\lambda$ variable). If the specifications for the Fourier basis are commented out, the information from the VMEC file is used.

### Running GVEC
To run the example, save the parameter file and pass it to the GVEC executable as:

```bash
  gvec_path/build/bin/gvec parameter.ini
```
  
It will generate for each `outputIter` a state file, the last one can then be used for the `convert_gvec_to_jorek` tool (also part of the GVEC repository).

## 3. Preparing for the JOREK simulation

### Running the GVEC to JOREK convertor

Before the initial conditions can be calculated, the GVEC output should be converted to a format that JOREK can read. This is done by using the `convert_gvec_to_jorek` utility:

```bash
gvec_path/build/bin/convert_gvec_to_jorek -r 41 -p 48 GVEC_W7A_UNSTABLE_1PA_mpol12_ntor10_r41_vmecBC_accGD_State_0000_00100000.dat gvec2jorek.dat
```

This utility constructs a non-axisymmetric flux-aligned grid with radial resolution given by the `-r` parameter and poloidal resolution given by the `-p` parameter, and then calculates finite element representations of the magnetic field and pressure from the GVEC equilibrium on this grid. The output file must always be named `gvec2jorek.dat`, as that is the file name which JOREK will try to open when looking for input from GVEC.

The `gvec2jorek.dat` file contains all the necessary information to reconstruct the equilibrium in JOREK except for the representation of the coil vacuum magnetic field, $\nabla \chi$. The vacuum field can be computed in one of two ways:
  - external to JOREK using a Dommaschk potential calculator.
  - internal to JOREK using the finite element basis and a Poisson solve. 

The next sections describe how to compute the Dommaschk potentials. If you do not intend to use these, you can skip to (4).

### Extracting the coil field with EXTENDER_P

The magnetic field representation used in the JOREK models requires that the part of the magnetic field that is generated by the coils be separated from the rest of the field, which is generated by plasma currents. Since a fixed-boundary equilibrium calculation provides only the total magnetic field, an intermediate step to extract the coil field is necessary. The [EXTENDER_P](https://dx.doi.org/10.1088/0029-5515/45/7/022) code will be used for this purpose.

The first step is to convert the output of GVEC to NetCDF format, which can then be read by EXTENDER_P. The `convert_gvec_to_vmec` utility, which is provided with the GVEC code can be used for this:
```bash
gvec_path/build/bin/convert_gvec_to_vmec -r 101 GVEC_W7A_UNSTABLE_1PA_mpol12_ntor10_r41_vmecBC_accGD_State_0000_00100000.dat wout_gvec_w7a_unst_1Pa.nc
```
where the `-r` parameter specifies the number of flux surfaces to include in the NetCDF file, the first parameter after that is the name of the GVEC output state file, and the last parameter is the desired name for the NetCDF file that will be created.

EXTENDER_P will output the $R$, $z$ and $\phi$ components of the coil field on an evenly spaced $(R,z,\phi)$ grid in MGRID format. Before running the code, one should create a text file called `suin` with the following contents:
```
nr 101
nz 101
nphi 91
```
This file specifies the resolution of the output grid. EXTENDER_P can then be run as follows:
```bash
EXTENDER_P -vmec_nyquist wout_gvec_w7a_unst_1Pa.nc -i suin -NU 150 -NV 150 -NI 7 -vacuumfield -write-mgrid -s gvec_w7a_unst_1Pa
```
Here, the `-vmec_nyquist` option specifies the name of the NetCDF file from the previous step, `-NU` and `-NV` give the resolution of the boundary surface on which the mirror currents are to be calculated, and `-NI` is the number of Newton-Cotes integration intervals. Finally, `-vacuumfield` specifies that the magnetic field due to the coils (the vacuum field, i.e. the field that would exist if the plasma was replaced by a vacuum) is to be calculated, `-write-mgrid` states that the output is to be written in MGRID format, and `-s` gives the suffix for the output file (in this case, the output file will be called `extended_field.gvec_w7a_unst_1Pa.nc`).

Note that EXTENDER_P is a pure MPI code, so there will be no benefit to assigning more than one core to each MPI task.

### Calculating the Dommaschk potential coefficients

Once EXTENDER_P finishes, the actual Dommaschk potential coefficients can be calculated using the `compute_dcoef_lsq.py` script, which is included the JOREK repository. The script can be run as follows:
```bash
jorek_path/util/compute_dcoef_lsq.py -R 1.99 -m 10 -l 9 -p 5 -o coef_1Pa_R1p99_l9.npy extended_field.gvec_w7a_unst_1Pa.nc
```
The script will do a least squares fit of the vacuum magnetic field from EXTENDER_P to the Dommaschk basis. One should specify the value of $R_0$ in the Dommaschk potential, which is the toroidally averaged $R$-coordinate of the vacuum field magnetic axis. However, instead of finding the vacuum magnetic axis, it is usually easier to use the $n=0$ Fourier mode of the total field axis (first number in `raxis`) as an initial guess, and then try to minimize the error from there. The `-m` parameter gives the number of toroidal modes, and, due to how the Dommaschk potentials are implemented in JOREK, must match the parameter `ntor` in the VMEC input namelist. The `-l` parameter gives the number of poloidal modes, and should be close to `mpol`, but does not need to match it exactly. The `-p` parameter is the periodicity of the device, and must match the `NFP` parameter in VMEC. The `-o` parameter specifies the name of the file to which the Dommaschk potential coefficients will be saved, and the last argument is the output file of EXTENDER_P in the previous step. More details about this calculation are available [here](assets/stellarator_setup/dommaschk_computation.pdf).

In order to compare the magnetic field given by the Dommaschk potential to the output of EXTENDER_P, the `compfld.py` script can be used:
```bash
jorek_path/util/compfld.py -R 1.9905208396958163 -F 4.720778214636723 -p 5 -o compare_1Pa_r50_l9.nc coef_1Pa_r50_l9.npy extended_field.gvec_w7a_unst_1Pa.nc
```
Note that this script depends on the `field.py`, `boundary.py` and `progress.py` modules from the [Zoidberg project](https://github.com/boutproject/zoidberg/tree/master/zoidberg). In `field.py`, the import statements for `boundary` and `update_progress` should be replaced with the following:
```python
import boundary
from progress import update_progress
```

The `compfld.py` script will calculate the volume-averaged square of the absolute error, $\langle |\vec B_2 - \vec B_1|^2\rangle$, as well as the volume-averaged squares of the errors relative to the EXTENDER_P field $\vec B_1$ and the Dommaschk potential field $\vec B_2$, $\langle |\vec B_2 - \vec B_1|^2/|\vec B_{1,2}|^2\rangle$. The `-R` and `-F` parameters must match the values of $R_0$ and $F_0$, respectively, as printed out by the `compute_dcoef_lsq.py` script, while `-p` is the periodicity. Finally, using the `-o` option with the desired file name allows one to create a NetCDF file which contains both the EXTENDER_P and Dommaschk fields for visual comparison in Paraview or Visit. The last two arguments are the names of the files containing the Dommaschk potential coefficients and the EXTENDER_P output, which were created in the previous steps.

The final step for the Dommaschk potential coefficients is to save them in a format that can be read by JOREK. This is done by saving the coefficients in a special namelist file, separate from the standard JOREK input file. The coefficients are held in the `dcoef` Fortran array, however the indices are ordered differently than in the NumPy array, and the Fortran array also skips unused coefficients (those whose toroidal mode number is not a multiple of the device periodicity). The following simple script can be used to write out the contents of the namelist file:
```python
import numpy as np

Nt = 10
Mp = 9
P = 5
coef = np.load("coef_1Pa_r50_l9.npy")
for i in range(Nt+1):
  for l in range(Mp+1):
    for j in range(4):
      if (coef[P*i,l,j] != 0.0):
        print(" dcoef(" + str(j+1) + "," + str(l) + "," + str(i) + ") = " + str(coef[P*i,l,j]))
```
The variables `Nt` and `Mp` should be set to the values passed to the `-m` and `-l` parameters, respectively, of the `compute_dcoef_lsq.py` script. The `P` variable is just the periodicity. The output of this script should then be copied to the namelist file. In addition, the namelist should include the `R_domm` variable, which should be set to the value of $R_0$, as printed out by `compute_dcoef_lsq.py`. An example of how the Dommaschk potential namelist file should look like is available here: [dc_W7A_unst_1Pa](assets/stellarator_setup/dc_W7A_unst_1Pa).

## 4. Reconstructing the equilibrium in JOREK using model 180

The next step is to compile JOREK with model 180, which is a "fake" model that does not evolve the plasma in time, but only calculates the initial conditions. The hard-coded options are to be set as follows: both `n_period` and `n_coord_period` should match `NFP` in the VMEC namelist, `n_coord_tor` should be `2*ntor + 1`, where `ntor` is the value of the corresponding parameter passed to VMEC.

If Dommaschk coefficients are used, `l_pol_domm` should be set to the number passed to the `-l` parameter in `compute_dcoef_lsq.py` in the previous step. If the code is run without Dommaschk potentials, `USE_DOMM=0` in jorek/defaults.mk and `l_pol_domm=0` must be set. This will mean that a Poisson solve is carried out as part of model 180 to compute $\chi$. Note that while this can take some time, it will make further calculations significantly faster, as computation of Dommaschk potentials is computationally demanding.

If 32-bit integers are used, an integer overflow of the matrix array indices may be reached for more complex setups due to the many DOF in the matrix solver. How switch to 64-bit integers is described here: [Compiling JOREK and solvers using long integers](../compiling/compilation_and_libraries/compiling_int64.md)

For the W7-A example considered here, using Dommaschk potentials, the hard-coded parameters have the following values:

```
model180
n_tor = 5
n_period = 5
n_coord_tor = 21
n_coord_period = 5
l_pol_domm = 9
n_plane = 40
n_vertex_max = 4
n_nodes_max = 15001
n_elements_max = 15001
n_boundary_max = 1001
n_pieces_max = 10001
```

After compilation, JOREK will need to be run with an input namelist. A sample namelist for the 1Pa example is provided below:
```fortran
&in1
 restart = .f.
 regrid  = .f.
 refinement=.false.
 tstep   = 1.
 nstep   = 1
 nout = 1

 freeboundary = .f.
 gvec_grid_import = .true.
 keep_current_prof = .f.
 fix_axis_nodes = .t.
 bc_natural_open = .t.

 F0    = 4.720778214636723

 rho_0 = 1.0
 rho_1 = 1.0
 rho_coef(1)  =  0.0
 rho_coef(2)  =  0.0
 rho_coef(3)  =  0.0
 rho_coef(4)  =  1.d0
 rho_coef(5)  =  5.d0
 T_0 = 1.d-6

 n_flux   = 41
 n_tht    = 48

 m_pol_bc = 24

 domm_file = "dc_W7A_unst_1Pa"

 iter_precon = 22
 gmres_m = 20
 gmres_4 = 1.d0
 gmres_max_iter = 200
 gmres_tol = 1.d-8
/
```

Here, `nstep` must always be 1, since there is no actual time stepping, and the initial conditions only need to be calculated once, while the value of `tstep` is irrelevant and does not affect the results. The `gvec_grid_import` parameter tells JOREK to look for the `gvec2jorek.dat` file, which must be present in the working directory, and import from there, while `bc_natural_open` allows the current to be nonzero on the boundary. Moreover, `F0` should be set to the value printed out by `compute_dcoef_lsq.py`. If run without Dommaschk potentials, `F0` is given by the `rbtor` variable in the NetCDF file produced by VMEC.

IMPORTANT: When running without the Dommaschk potentials, a density profile must also be provided. Otherwise, the density and the temperature will return NaN values. 

`n_flux` and `n_tht` must match the values passed to the `-r` and `-p` parameters of `convert_gvec_to_jorek`, respectively, and `domm_file` gives the path to the namelist where the Dommaschk potential coefficients from the previous step are stored. The file used in this tutorial is available here: [dc_W7A_unst_1Pa](assets/stellarator_setup/dc_W7A_unst_1Pa).

JOREK can then be run as follows:
```bash
mpirun -n 3 ./jorek_model180 < intear
```
where the number passed to the `-n` parameter must be a multiple of `(n_tor+1)/2`, and `intear` is the name of the namelist file given above. This run will produce, among others, an output file called `jorek_restart.h5`, which contains the initial conditions for the reduced MHD variables. This file will be required in the next step.

## 5. Running the simulation in JOREK in time using model 183

The actual time stepping is done in model 183. Keeping the same hard-coded parameters as before, JOREK should be recompiled with model 183. In the example considered here, a tearing mode will be simulated. The first part of the simulation will proceed with five-fold periodicity to allow the plasma to re-equilibrate as small residual forces are damped out. The following namelist file can be provided to JOREK:
```fortran
&in1
 restart = .t.
 regrid  = .f.
 refinement=.false.
 tstep_n = 1. 10. 100.
 nstep_n = 20 20  5
 nout = 1

 freeboundary = .f.
 keep_current_prof = .t.
 init_current_prof = .t.
 time_evol_scheme = "implicit Euler"
 fix_axis_nodes = .t.
 bc_natural_open = .f.
 eta_T_dependent = .f.
 visco_T_dependent = .f.
 zkpar_T_dependent = .f.

 F0    = 4.720778214636723

 rho_0 = 1.0
 rho_1 = 1.0
 rho_coef(1)  =  0.0
 rho_coef(2)  =  0.0
 rho_coef(3)  =  0.0
 rho_coef(4)  =  1.d0
 rho_coef(5)  =  5.d0
 T_0 = 1.256637061435917E-006

 n_flux   = 41
 n_tht    = 48

 eta   = 1.d-6
 visco = 0.d-9

 eta_num = 0.d-19
 visco_num = 1.d-15

 D_par  = 0.d0
 D_perp = 0.d-7

 ZK_par  = 0.d0
 ZK_perp = 0.d-7

 domm_file = "dc_W7A_unst_1Pa"

 heatsource     = 0.d0
 particlesource = 0.d0

 iter_precon = 22
 gmres_m = 1000
 gmres_4 = 1.d0
 gmres_max_iter = 400
 gmres_tol = 1.d-8
/
```
Here, `restart` is true, specifying that the simulation is to be initialized from the `jorek_restart.h5` file  produced by the previous JOREK run (initial conditions step). In addition, `init_current_prof` must be set to true, so that the current source is set to the initial condition for the current. Finally, the implicit Euler method should be used for temporal discretization (`time_evol_scheme`), as it is more dissipative than the Crank-Nicolson method, and will damp out oscillations due to residual forces relatively quickly. The Gears method can also be used.

JOREK can be run as follows:
```bash
mpirun -n 3 ./jorek_model183 < intear
```
where, as before, the number passed to the `-n` parameter must be a multiple of `(n_tor+1)/2`, and `intear` is the name of the namelist file given above. Again, this run will produce a file called `jorek_restart.h5`, which will be used to initialize the next JOREK run.

The final step is to expand the simulation domain to the full torus to allow a tearing mode to develop. This requires setting the hard-coded parameters `n_tor` to 21 and `n_period` to 1, and then recompiling JOREK. After recompilation, the following namelist file can be used for this final step:
```fortran
&in1
 restart = .t.
 regrid  = .f.
 refinement=.false.
 tstep   = 25.
 nstep   = 160
 nout = 1

 freeboundary = .f.
 keep_current_prof = .t.
 init_current_prof = .f.
 fix_axis_nodes = .t.
 bc_natural_open = .f.
 eta_T_dependent = .f.
 visco_T_dependent = .f.
 zkpar_T_dependent = .f.

 F0    = 4.720778214636723

 rho_0 = 1.0
 rho_1 = 1.0
 rho_coef(1)  =  0.0
 rho_coef(2)  =  0.0
 rho_coef(3)  =  0.0
 rho_coef(4)  =  1.d0
 rho_coef(5)  =  5.d0
 T_0 = 1.256637061435917E-006

 n_flux   = 41
 n_tht    = 48

 eta   = 1.d-6
 visco = 0.d-9

 eta_num = 0.d-19
 visco_num = 1.d-15

 D_par  = 0.d0
 D_perp = 0.d-7

 ZK_par  = 0.d0
 ZK_perp = 0.d-7

 domm_file = "dc_W7A_unst_1Pa"

 heatsource     = 0.d0
 particlesource = 0.d0

 iter_precon = 22
 gmres_m = 1000
 gmres_4 = 1.d0
 gmres_max_iter = 400
 gmres_tol = 1.d-8
/
```
This file is almost the same as in the previous step, except for the step size (`tstep`) and number of steps (`nstep`), as well as `init_current_prof` being set to false (the current source will be set to the value of the current source stored in `jorek_restart.h5`, instead of the value of current) and `time_evol_scheme` not being specified (the Crank-Nicolson method is the default).

JOREK can then be run as follows:
```bash
mpirun -n 11 ./jorek_model183 < intear
```
where, as before, `intear` is the name of the namelist file given above.
