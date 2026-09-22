---
title: "Normalization"
nav_order: 2
parent: "Physics Models"
layout: default
render_with_liquid: false
---
# Normalization in JOREK

* Subscript "SI" denotes quantities in SI units while quantities in JOREK units are written without subscripts.
* Note, that the important normalization factors $\sqrt{\mu_{0}/\rho_{0}}$ and $\sqrt{\mu_{0}\rho_{0}}$ are written into the JOREK logfile for convenience. 
* Note, that you need to set the input parameters `central_density` and `central_mass` for these values to be correct. 

## Connection between SI and normalized units

| Quantity | Physical Unit (SI) | Connection | Description |
| :--- | :--- | :--- | :--- |
| Major radius | $R_{SI}[m]$ | $=R$ | Major radius |
| Vertical coordinate | $Z_{SI}[m]$ | $=Z$ | Vertical coordinate |
| Magnetic field vector | $B_{SI}[T]$ | $=B$ | Magnetic field vector |
| Electric field vector | $E_{SI}[Vm^{-1}]$ | $=E/\sqrt{\mu_{0}\rho_{0}}$ | Electric field vector |
| Poloidal magnetic flux | $\Psi_{SI}[Tm^{2}]$ | $=\Psi$ | Poloidal magnetic flux |
| Toroidal current density | $j_{\phi,SI}[Am^{-2}]$ | $=-j/(R~\mu_{0})$ | Toroidal current density; $j_{\phi,SI}=j_{SI}\cdot\hat{e}_{\phi}$ |
| Runaway electron number density | $n_{r,SI}[m^{-3}]$ | $=n_{r}(\frac{1}{eR})\sqrt{\frac{\rho_{0}}{\mu_{0}}}$ | Runaway electron number density |
| Runaway electron parallel momentum | $P_{\parallel,SI}[kg~m~s^{-1}]$ | $=P_{\parallel}m_{e0}c$ | Runaway electron parallel momentum |
| Particle density | $n_{SI}[m^{-3}]$ | $=\rho~n_{0}$ | Particle density ($\rho$ is the normalized density profile, which should be given as an input) |
| Impurity number density | $n_{imp,SI}[m^{-3}]$ | $=\rho_{imp}\;n_{0}\;\mu_{imp}$ | Impurity number density |
| Mass density | $\rho_{SI}[kg~m^{-3}]$ | $=\rho~\rho_{0}$ | Mass density = ion mass X particle density |
| Impurity mass density | $\rho_{imp,SI}[kg/m^{3}]$ | $=\rho_{imp}\rho_{0}$ | Impurity mass density |
| Temperature | $T_{SI}[K]$ | $=T/(k_{B}\mu_{0}n_{0})$ | Temperature electron + ion temperature |
| Temperature (eV) | $T_{eV}[eV]$ | $=T/(e~\mu_{0}n_{0})$ | Temperature in eV |
| Poloidal current stream function | $FF_{SI}^{\prime}[T^{2}m^{2}/(Weber/rad)]$ | $=FF^{\prime}$ | Poloidal current stream function $F=RB_{\phi}$ ${}^{\prime}=d/d\psi$ |
| Plasma pressure | $p_{SI}[Nm^{-2}]$ | $=\rho~T/\mu_{0}$ | Plasma pressure |
| Velocity vector | $v_{SI}[ms^{-1}]$ | $=v/\sqrt{\mu_{0}\rho_{0}}$ | Velocity vector |
| Parallel velocity component | $v_{\parallel,SI}[ms^{-1}]$ | $=v_{\parallel}\cdot B_{SI}/\sqrt{\mu_{0}\rho_{0}}$ | Parallel velocity component, where $B_{SI}=\vert B_{SI}\vert $ |
| Velocity stream function | $u_{SI}[ms^{-1}]$ | $=u/\sqrt{\mu_{0}\rho_{0}}$ | $Ru$ is the velocity stream function, $F_0 u$ is the potential |
| Toroidal vorticity | $\omega_{\phi,SI}[m^{-1}s^{-1}]$ | $=\omega/\sqrt{\mu_{0}\rho_{0}}$ | Toroidal vorticity |
| Time | $t_{SI}[s]$ | $=t\cdot\sqrt{\mu_{0}\rho_{0}}$ | Time |
| Growth rate | $\gamma_{SI}[s^{-1}]$ | $=\gamma/\sqrt{\mu_{0}\rho_{0}}$ | Growth rate; $\gamma_{SI}=\ln[E_{SI}(t_{2})/E_{SI}(t_{1})]/[2\Delta t_{SI}]$ Energy $E_{SI}[J]$ |
| Resistivity | $\eta_{SI}[\Omega m]$ | $=\eta\cdot\sqrt{\mu_{0}/\rho_{0}}$ | Resistivity, see also notes on Spitzer resistivity |
| Hyper-resistivity | $\eta_{num,SI}[\Omega m^{2}]$ | $=\eta_{num}\cdot\sqrt{\mu_{0}/\rho_{0}}$ | Hyper-resistivity |
| Dynamic viscosity | $\mu_{SI} [kg~m^{-1}s^{-1}]$ | $=\mu\cdot\sqrt{\rho_{0}/\mu_{0}}$ | Dynamic viscosity |
| Hyper-viscosity | $\mu_{num,SI}[kg~ms^{-1}]$ | $=\mu_{num}\cdot\sqrt{\rho_{0}/\mu_{0}}$ | Hyper-viscosity |
| Kinematic viscosity | $\nu_{SI}[m^{2}s^{-1}]$ | $=\mu_{SI}/\rho_{SI}$ | Kinematic viscosity ($\rho_{SI}$ is the local mass density in $kg~m^{-3}$) |
| Particle diffusivity | $D_{SI}[m^{2}s^{-1}]$ | $=D/\sqrt{\mu_{0}\rho_{0}}$ | Particle diffusivity (\parallel or $\perp$); Usually, $D_{\parallel}=0$ |
| Heat diffusivity | $K_{SI}[kg~m^{-1}s^{-1}]$ | $=K\cdot\sqrt{\rho_{0}/\mu_{0}}/(\gamma-1)$ | Heat diffusivity (\parallel or $\perp$), where $\chi_{SI} [m^{2}s^{-1}]=K_{SI}/\rho_{SI}$ and $\kappa_{SI} [m^{-1}s^{-1}]=n_{SI}\chi_{SI}$ |
| Heat source | $S_{T,SI}[Wm^{-3}]$ | $=S_{T}/((\gamma-1)\mu_{0}\sqrt{\mu_{0}\rho_{0}})$ | Heat source |
| Particle source | $S_{\rho,SI}[kg~s^{-1}m^{-3}]$ | $=S_{\rho}\cdot\sqrt{\rho_{0}/\mu_{0}}$ | Particle source |
| Wall resistivity | $\eta_{wall,thin,SI} [\Omega]$ | $=\eta_{wall,thin}\cdot\sqrt{\mu_{0}/\rho_{0}}$ | Wall resistivity (relevant for JOREK-STARWALL); $\eta_{wall,thin,SI} [\Omega] = \eta_{wall,SI} [\Omega m] / d_{wall} [m]$. Example ITER: $8\cdot10^{-7}\Omega m / (6cm) = 1.33\cdot10^{-5}\Omega$ |
| Ionisation/recombination rate | $R_{ion/rec,SI}[m^{-3}s^{-1}]$ | $=R_{ion/rec}/(\sqrt{\mu_{0}\rho_{0}}n_{0})$ | Ionisation and recombination rate |
| Ionisation energy | $E_{ion,SI}[J]$ | $=\xi_{ion}/((\gamma-1)\mu_{0}n_{0})$ | Ionisation energy |
| Radiation rate | $L_{rad,SI}[Wm^{3}]$ | $=L_{rad}/((\gamma-1)\mu_{0}\sqrt{\mu_{0}\rho_{0}}n_{0}^{2}\frac{m_{i}}{m_{imp}})$ | Radiation rate (model501) |
| Radiation power density | $P_{rad,SI}[Wm^{-3}]$ | $=P_{rad}/((\gamma-1)\mu_{0}\sqrt{\mu_{0}\rho_{0}})$ | Radiation power density (model501) |
| Particle charge | $q_{SI}[As]$ | $=q\sqrt{\rho_{0}/\mu_{0}}$ | Particle charge |
| Neoclassical friction rate | $\mu_{neo,SI}[s^{-1}]$ | $=\mu_{neo}/\sqrt{\rho_{0}\mu_{0}}$ | Neoclassical friction rate |

## Having defined:

* **$n_{0}[m^{-3}]$**: $= \text{central\_density } 10^{20}$. central_density gets a default value in `preset_parameters.f90` and should be specified in the input file.
* **$\rho_{0}[kg~m^{-3}]$**: $= \text{central\_mass}*m_{AMU}*n_0$. central_mass gets a default value in `preset_parameters.f90` and should be specified in the input file. 
* **$\gamma$**: $= \text{GAMMA}$. GAMMA gets a default value of $5/3$ in `preset_parameters.f90`. 
* **$\mu_{imp}$**: `m_i_over_m_imp` $= m_{i}/m_{imp}$. Defined by the impurity species.

## Useful Constants:

* $m_{e0}=0.911\cdot10^{-30}kg$
* $m_{AMU}=1.661\cdot10^{-27}kg$
* $m_{deuterium}=2.014101777811AMU$
* $m_{tritium}=5.007\cdot10^{-27}kg$
* $\mu_{0}=4\cdot\pi\cdot10^{-7}Vs/(Am)$
* $e=1.602176565\cdot10^{-19}C$
