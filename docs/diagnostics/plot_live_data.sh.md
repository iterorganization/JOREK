---
title: "plot_live_data.sh script"
nav_order: 16
render_with_liquid: false
parent: "Running JOREK"
---

- Run the script from your simulation folder to plot for instance the **energy time traces** during or after a simulation.
- For example (click to enlarge): ![](energies.png?400)
- The script requires a recent gnuplot version (on many machines available directly or via the 'module' command)
- It is also possible to plot **growth rates** (`plot_live_data.sh -q g`), simulation times associated with the **time steps** (`plot_live_data.sh -q t`), and **input profiles** (`plot_live_data.sh -q i`) and can easily be extended in the future
- **List all available quantities with `plot_live_data.sh -l`**

## Usage Information (plot_live_data.sh -h)

  ```
  Usage: plot_live_data.sh [-h -f <file> -l -q <qtty> -ps -(no)log]

    -h         Print this usage information
    -f <file>  Take data from <file> instead of 'macroscopic_vars.dat'
    -l         List plottable quantities
    -q <qtty>  Plot the given quantity (default: '-q energies')
    -ps        Plot to .ps files (default: plot to screen)
    -(no)log   (Non-)Logarithmic y-axis (default: '-log')

  Remarks:
  - The beginning of a quantity name is sufficient, e.g., the command
      'plot_live_data.sh -q gr' will plot growth_rates.

  ```
## Check energy conservation

### Total energy conservation

Live data allows you to check how JOREK conserves energy (or not) in real time. **The energy conservation diagnostic must be adapted to each model** since the addition of new physical terms can lead to new energy channels. Considering a basic visco-resistive MHD model, the evolution of the total energy density is

\begin{equation}
\label{eq:energy_density}
\partial_t w + \nabla\cdot\mathbf{\Gamma}=\tau_{nc} 
\end{equation}

where the energy density and the energy flux (omitting $\mu_0$ factors) are respectively

\begin{align}
w &\equiv \frac{\rho}{2}|\mathbf{v}|^2 + \frac{p}{\gamma-1}+\frac{|\mathbf{B}|^2}{2}   
 \mathbf{\Gamma} &\equiv \left[ \frac{\rho}{2}|\mathbf{v}|^2  + \frac{\gamma}{\gamma-1}p \right]\mathbf{v} + \frac{\mathbf{q}}{\gamma-1} + \mathbf{E}\times\mathbf{B} 
\end{align}

the RHS of Eq. \eqref{eq:energy_density} ($\tau_{nc}$) represents non-conservative terms such as sources (such as heating sources) and sinks (such as radiation). Dissipative terms such as resistive and viscous terms should in principle not appear inside ($\tau_{nc}$). However if for example, $\eta \neq \eta_{ohmic}$, some of the loss magnetic energy will not be converted into thermal energy through the Ohmic heating term, leading to a non-conservative term. The derivation of \eqref{eq:energy_density} is sketched in the following ![slides ](artola_energy_conservation_diagnostics.pdf). The energy conservation diagnostic is based on the volume integral of Eq. \eqref{eq:energy_density} over the JOREK domain, which is

\begin{equation}
\label{eq:total_energy_ev}
\partial_t E_{tot} = \textrm{Boundary fluxes} + \textrm{Non-conservative/source terms}  
\end{equation}

and

\begin{align}
E_{tot} &\equiv \int_V w dV   
 \textrm{Boundary fluxes} &\equiv -\oint_S \mathbf{\Gamma}\cdot\mathbf{n}dS  
\textrm{Non-conservative/source terms} &\equiv \int_V \tau_{nc} dV
\end{align}

You can check that Eq. \eqref{eq:total_energy_ev} is satisfied by checking the LFS and RHS terms with the command

  jorek_folder/util/plot_live_data.sh -q energy_conservation

In case of good energy conservation, both lines should overlap as seen in the example below.

![](example_energy_conservation1.png?nolink)

You can also plot the different energy types and their time derivatives with the commands

  ```
  jorek_folder/util/plot_live_data.sh -q integrated_energies
  jorek_folder/util/plot_live_data.sh -q dEdt

  ```
which are named with the following labels

\begin{align}
\textrm{Total energy} &\equiv E_{tot}   
\textrm{Magnetic} &\equiv \int_V \frac{|\mathbf{B}|^2}{2} dV   
\textrm{Thermal energy} &\equiv \int_V \frac{p}{\gamma-1} dV   
\textrm{Kinetic parallel} &\equiv \int_V\frac{\rho}{2}|\mathbf{v}_\parallel|^2 dV   
\textrm{Kinetic perpendicular} &\equiv \int_V\frac{\rho}{2}|\mathbf{v}_u|^2 dV 
\end{align}
where $\mathbf{v}_u=-R^2\nabla u\times \nabla\phi$.

You can also plot different boundary fluxes with the command

  jorek_folder/util/plot_live_data.sh -q bnd_fluxes

which are labelled as

\begin{align}
\textrm{"p vn"} &\equiv \int_S \frac{\gamma}{\gamma-1}p \mathbf{v}\cdot\mathbf{n}dS   
\textrm{"kinpar-flux"} &\equiv \int_S \frac{\rho}{2}|\mathbf{v}|^2 \mathbf{v}\cdot\mathbf{n}dS    
\textrm{"qn-par"} &\equiv  \int_S  \frac{\mathbf{q}_{\parallel}\cdot\mathbf{n}}{\gamma-1}dS   
\textrm{"qn-perp"} &\equiv  \int_S  \frac{\mathbf{q}_{\perp}\cdot\mathbf{n}}{\gamma-1}dS   
\end{align}

Note that if you set sheath boundary conditions in all the boundary of the domain (**not just the divertor**), the sum of the 4 previous fluxes corresponds to the energy flowing through the sheath and the formula provided by Stangeby (see [Running with sheath boundary conditions for the heat-flux](/JOREK/howto/physics_options/boundary_conditions/sheath_heatflux_bc.html)).

The Poynting flux is typically 0 in fixed boundary simulations (see next section on magnetic energy conservation how to plot it.)

In the current JOREK version (adapted for RMHD) the following terms are considered in the diagnostic as non-conservative terms

\begin{equation}
\textrm{Non-conservative/source terms} = \int_V \left( \frac{S_T}{\gamma-1} + S_{mag} + (\eta_{ohmic}-\eta) |\mathbf{J}|^2  -\boldsymbol{\nu}_\parallel\cdot\mathbf{v}_\parallel  \right) dV
\end{equation}

where $S_T = S_{T_i}+S_{T_e}$ is the heating source for ions and electrons together (used in the temperature equations), $S_{mag}$ is the magnetic energy source provided by the current source term (see next section) and the last term is the energy dissipated through parallel viscosity. These terms can be plotted with

  jorek_folder/util/plot_live_data.sh -q dissipative_terms

Note that all of these quantities can be **recovered via post-processing** with [jorek2_postproc](jorek2_postproc.md) with the command "expressions_int". See usage with "help expressions_int" and list all the available quantities (about 67) by typing "expressions_int" in jorek2_postproc.



### Magnetic energy conservation

Similar to the total energy equation, the different plots of the magnetic energy equation can be plotted with the command

  ```
  jorek_folder/util/plot_live_data.sh -q mag_energy_balance

  ```
![](screenshot_from_2020-04-19_13-20-17.png?800)

The equation for the magnetic energy and the plotted terms can be found in the document ![](mag_energy_evolution.pdf)  (also you can see the corresponding ![ pull request](https://git.iter.org/projects/STAB/repos/jorek/pull-requests/376/overview) )
