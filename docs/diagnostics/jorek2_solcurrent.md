---
title: "jorek2_SOLcurrent"
nav_order: 31
render_with_liquid: false
parent: "General"
---

## jorek2_SOLcurrent

This diagnostics postprocesses jorek_restart.h5 to determine scrape-off layer currents (thermo electric current between target plates) by solving Stangeby eq 17.29. Run it as 

```
jorek2_SOLcurrent < input
```

with "input" the input namelist file used for running jorek (jorek_model600 < input). The output it generates will be put into a folder "SOLcurrent" (which will be generated if it does not exist yet), and the file is called "stepXXXXXX.dat" (e.g. step000200.dat if jorek_restart.h5 was jorek00200.h5). The content of the output file is the solution of the equation per fluxline along the boundary with some intermediate output. Example output (truncated to two output lines rather than the n_bnd_elements output lines in reality):

```
 #            R_b               Z_b             Phi_b               R_e               Z_e             Phi_e              ne_b              ne_e              Te_b              Te_e                 L      av_sigma_par         P_h - P_c           pe_term         j_hat_par             j_par             bdotn            j_wall                dl             I_loc
    5.56420000E+00   -4.55671484E+00    0.00000000E+00    5.56420000E+00   -4.55671484E+00    0.00000000E+00    5.82638640E+16    5.82638640E+16    4.51428697E-01    4.51428697E-01    0.00000000E+00    0.00000000E+00    0.00000000E+00    0.00000000E+00    0.00000000E+00    0.00000000E+00   -5.36382742E-02    0.00000000E+00    1.48283537E-02    0.00000000E+00
    5.56420000E+00   -4.54188649E+00    0.00000000E+00    4.16796881E+00   -3.90608038E+00   -6.46920674E+00    4.45768773E+17    4.09077666E+17    4.51369085E-01    5.35422967E-01    3.17676936E+01    3.72828298E+03   -8.04289104E-02    3.90718320E-01    8.92291495E-02    4.20341276E+01   -5.26322395E-02   -1.93477645E+00    1.51329629E-02   -1.02361670E+00
```

_b means at the beginning of the fieldline for which the equation is solved, _e at the end of the fieldline, R,Z and Phi are the coordinates, ne, Te electron density and temperature, L connection length, av_sigma_par the average parallel conductivity, P_h-P_c the pressure drop along the fluxtube from the hot to cold end (this can be used as sanity check by comparing the pressures at the (R_b,Z_b,Phi_b) location to the (R_e,Z_e,Phi_e) location), pe_term the pressure integral term in the equation (the last term), j_hat_par the solution to the equation, j_par the physical thermo-electric parallel current density at location (R_b,Z_b,Phi_b), bdotn the pitch angle dot product at (R_b,Z_b,Phi_b), j_wall the thermo-electric current density on the wall (so accounting for pitch angle), dl the poloidal extent of this flux tube, I_loc the current onto the wall of this flux tube.

Running the diagnostics takes about ~3 minutes on a single omp core if jorek2_SOLcurrent is compiled without DEBUG options for a production FE grid of n_bnd_elements ~ 500, or less if using omp parallelisation (set "export OMP_NUM_THREADS=..." before calling the diagnostic). While running the diagnostic the output will let you know where in the process it is, usually the "Initialising element neighbours" takes about half the time, with the other half being the actual fieldline tracing and solving for the equation (of which you will receive intermediate progress updates as "starting point =         49 of        482")

When the diagnostic is finished, it will let you know the sum of all I_loc values as a sanity check:

```
Sanity check: total current over the whole wall (i.e. discretisation error, usually of order 10 A, should be negligible compared to target current, usually of order 10 kA) =    -0.11341611E+02 A
```

In principle this sum should go to 0 as the amount of current per flux tube should be the same at the inner and outer target, but with opposite sign.


To calculate the actual target currents, use any programming language of your choice to postprocess the stepXXXXXX.dat file based on what extent of the currents you want to take into account. Example code in python below:

```python
#importing modules
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline

#load .dat file
filepath = '/path/to/sim/folder/.../SOLcurrent/stepXXXXXX.dat'
j_SOL_dat = np.genfromtxt(filepath,dtype=None).T

#definitions of all variables
yl_SOLcurrentdat = ["$R_b$ (m)","$Z_b$ (m)","$\phi_b$","$R_e$ (m)","$Z_e$ (m)","$\phi_e$","$n_{e,b}$ (m⁻³)",
              "$n_{e,e}$ (m⁻³)","$T_{e,b}$ (eV)","$T_{e,e}$ (eV)","$L$ (m)",
              "$\sigma_{\parallel,av}$ ($\Omega$⁻¹m⁻¹)","$P_{e,h} - P_{e,c}$ (Pa)","$p_e$ term",
              "$\hat{j}_{\parallel}$ (-)","$j_{\parallel}$ (A/m²)", "$\hat{b} \cdot \hat{n}$", 
              "$j_{wall}$ (A/m²)","dl (m)", "$I_{loc}$ (A)"]

#target mask (ITER)
[R,Z] = j_SOL_dat[0:2]
[R2,Z2] = j_SOL_dat[3:5]
mask_IT = (R < 5) & (Z < -2.5) & (-0.625*R + Z > -6.5)
mask_OT = (R > 5.55) & (Z < -3.3) & (Z > -4.5)

#SOL current at both targets
I_SOL_IT = sum(j_SOL_dat[-1][mask_IT])
I_SOL_OT = sum(j_SOL_dat[-1][mask_OT])

print(f'sum I_loc = {(sum(j_SOL_dat[-1])/1e3):.2f} kA')
print(f'I_SOL into   IT = {(-I_SOL_IT)/1e3:.1f} kA')
print(f'I_SOL out of OT = {(I_SOL_OT) /1e3:.1f} kA')

#plot quantities
plot_from = 10
for j in range(len(j_SOL_dat[plot_from:])):
    i = j + plot_from
    qty = j_SOL_dat[i]
    plt.figure()
    plt.plot(Z[mask_IT],qty[mask_IT],'.',label='started at IT')
    plt.plot(Z2[mask_OT],qty[mask_OT],'x',label='ended at IT')
    plt.xlabel('Z (IT end of fieldline)')
    plt.ylabel(yl_SOLcurrentdat[i])
    plt.xlim([-4,-3])
    plt.grid("show")
    plt.legend()
    plt.show()
```
This prints the currents at the two targets and the sum over all currents (the sanity check mentioned earlier):
```
sum I_loc = -0.01 kA
I_SOL into   IT = 15.5 kA
I_SOL out of OT = 15.5 kA
```
And generates a number of plots, such as:

![](j_par_jorek2_solcurrent.png?400)

In this figure it can be seen that the calculation for tracers starting at the inner target are the same as for tracers starting at the outer target, as should be the case.

On the sign of j in this routine: current is defined as positive if it comes out of the wall at the starting position (R_b,Z_b,Phi_b). Usually the outer target is hotter, so current usually flows from OT to IT, giving positive sign at the OT and negative sign at the IT.

A simple bash script could be used to loop over many restart files, to avoid the labour of running it manually for each restart file separately.

For more information on implementation, results on an ITER case and comparison to SOLPS results of SOL current, see [Daniël's Msc. Thesis report](https://research.tue.nl/files/340298012/1367803_-Maris_D.-MSc_thesis_Thesis-_NF.pdf) section 4.1
