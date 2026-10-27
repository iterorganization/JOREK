---
title: "Model 600 weak form"
nav_order: 3
parent: "Reduced MHD models"
grand_parent: "Base Fluid Models"
layout: default
render_with_liquid: false
---

# Model 600: equations in weak form

This page documents the equations of the model-600 element routine
`models/model600/mod_elt_matrix_fft.f90` in weak form. It covers the equations
for $\psi$, $u$, $j$, $\omega$, $\rho$, $v_\parallel$, $\rho_{imp}$, $\rho_n$,
$T_i$, $T_e$ and $T$; see [Status](#status) for what is not covered yet.

The equations are written, as Python, in
`util/equation_codegen/src/jorek_equations/model600.py`, and the equation
checker compares them with the Fortran term by term. All math on this page is
generated from that file (see [Regenerating this page](#regenerating-this-page)).
Where the Python and the Fortran differ, the difference is known and recorded:
it is listed under [Known differences](#known-differences) and explained in
[`JOREK_FINDINGS.md`](https://github.com/iterorganization/JOREK/blob/develop/util/equation_codegen/JOREK_FINDINGS.md).

**The integral is implied.** Each equation is multiplied by a test function
$v$ (a finite-element basis function, not a velocity) and integrated over the
plasma volume, $\int \ldots\,\mathrm{d}V$ with
$\mathrm{d}V = w_V\,\mathrm{d}s\,\mathrm{d}t_{\mathrm{el}}\,\mathrm{d}\phi$ and
$w_V = R\,\mathcal{J}$, where $(s, t_{\mathrm{el}})$ are the element
coordinates and $\mathcal{J}$ their Jacobian. (The code calls the second
element coordinate `t`; on this page $t$ is always time.) The integral, the
common volume weight $w_V$ and the coordinate differentials are not written,
so $\frac{\partial}{\partial t}\big(v\,\rho\big) = v\,S_p + \ldots$ stands for
$\frac{\partial}{\partial t}\int v\,\rho\,\mathrm{d}V = \int \big(v\,S_p + \ldots\big)\,\mathrm{d}V$.
Terms with derivatives of $v$ come from integration by parts; see
[Conventions](#conventions) for the assumptions behind this.

$\nabla_{\mathrm{pol}}$ is the poloidal gradient and $[a,b]$ the poloidal
bracket. The term groups written $\mathcal{C}$, $\mathcal{D}$, $\mathcal{T}$,
$Q$, ... are defined under [Shared operators](#shared-operators), and the
other symbols under [Notation](#notation). Refer also to the
[RMHD model](rmhd_model.md), [notation](../../notation.md) and
[normalization](../../normalization.md) pages.

## Variables

The model is configured at compile time by the logical parameters in
`models/model600/mod_model_settings.f90` (set with `util/config.sh`). When a
variable is switched off, the element routine replaces it by a constant:
$\rho = 1$, $v_\parallel = 0$, $\rho_n = 0$, $\rho_{imp} = 0$, and its
equation is not assembled. The equations below are written with every variable
present.

| Python | Symbol | Meaning | Present when |
|---|---|---|---|
| `psi` | $\psi$ | poloidal magnetic flux | always |
| `u` | $u$ | stream function of the $E\times B$ velocity | always |
| `j` | $j$ | toroidal current, $j = -R\,\mathbf{j}\cdot\mathbf{e}_\phi = \Delta^*\psi$; it includes a factor $R$ (see the [RMHD model](rmhd_model.md)) | always |
| `omega` | $\omega$ | vorticity, $\omega = \nabla^2_{\mathrm{pol}} u$ (the [vorticity definition](#eq-w) below) | always |
| `rho` | $\rho$ | total mass density, impurities included, in units of the main-ion mass | `with_rho = .true.` |
| `rhoimp` | $\rho_{imp}$ | impurity mass density, in the same units | `with_impurities = .true.` |
| `rhon` | $\rho_n$ | neutral density | `with_neutrals = .true.` |
| `vpar` | $v_\parallel$ | parallel velocity per unit field: the parallel flow is $v_\parallel\mathbf{B}$ | `with_vpar = .true.` |
| `T` | $T$ | total temperature $T_i + T_e$ (single-temperature model, where $T_i = T_e = T/2$) | `with_TiTe = .false.` |
| `Ti` | $T_i$ | ion temperature (two-temperature model) | `with_TiTe = .true.` |
| `Te` | $T_e$ | electron temperature (two-temperature model) | `with_TiTe = .true.` |

Some quantities take a different temperature in the two models
(`with_TiTe` selects the two-temperature model):

| Python | Symbol | Two-temperature model | Single-temperature model |
|---|---|---|---|
| `T_or_Te` | $\check{T}_e$ | $T_e$ | $T$ |
| `Te_gen` | $\hat{T}_e$ | $T_e$ | $T/2$ |
| `ion_temperature` | $\tilde{T}_i$ | $T_i$ | $T/2$ |

In the single-temperature model this distinction is deliberate: the rate and
charge-state functions ($S_{ion}$, $S_{rec}$, $\alpha_e$, $\eta$, ...) take
$\check{T}_e = T$, while the electron pressure of the induction equation and
the ion pressure of the diamagnetic terms take $T/2$. Do not "correct" one to
match the other.

Derived quantities that appear in several equations:

| Symbol | Definition | Meaning |
|---|---|---|
| $\rho_{main}$ | $\rho - \rho_{imp}$ | main-ion mass density |
| $n_i$ | $\rho + \alpha_i\,\rho_{imp}$ | ion density, $\alpha_i = m_i/m_{imp} - 1$ |
| $n_e$ | $\rho + \alpha_e(\check{T}_e)\,\rho_{imp}$ | electron density, $\alpha_e = (m_i/m_{imp})\,Z_{imp} - 1$ |
| $n_e^{\mathrm{c}}$ | $\rho^{\mathrm{c}} + \alpha_e(\check{T}_e)\,\rho_{imp}^{\mathrm{c}}$ | electron density from the corrected densities (not $\mathrm{corr}(n_e)$) |
| $x^{\mathrm{c}}$ | $\mathrm{corr}(x)$ | density $x$ corrected for negative values (`corr_neg_dens`); $x^{\mathrm{c}} = x$ where $x$ is well above zero |
| $p_e(\Theta)$ | $\rho\,\Theta + \rho_{imp}\,A_e(\Theta)$ | electron pressure, $A_e(\Theta) = \alpha_e(\Theta)\,\Theta$; $\Theta = T_e$ in the electron energy equation and $\Theta = \hat{T}_e$ in the induction equation |
| $p_i$ | $n_i\,\tilde{T}_i$ | ion pressure; written $p_i^{dia}$ in the diamagnetic terms |
| $p$ | two-temperature: $\rho\,(T_i + T_e) + \rho_{imp}\,\big(\alpha_i\,T_i + A_e(T_e)\big)$; single-temperature: $\rho\,T + \rho_{imp}\,A_{imp}(T)$ | total pressure, $A_{imp}(T) = \alpha_{imp}(T)\,T$ with $\alpha_{imp} = \frac{1}{2}(m_i/m_{imp})(Z_{imp}+1) - 1$ |
| $W_{ion}$ | $E_{ion}(\check{T}_e)\,\rho_{imp} + E_{ion}^{bg}\,(\rho - \rho_{imp})$ | ionization potential energy of the impurities and of the main ions |
| $B^2$ | $(F_0^2 + \lvert\nabla_{\mathrm{pol}}\psi\rvert^2)/R^2$ | square of the magnetic field |

## Equations

<!-- BEGIN GENERATED: equations -->

<a id="eq-psi"></a>

### Induction equation (`var_psi`)

$$
\begin{aligned}
\frac{\partial}{\partial t}\Big(\frac{1}{R^{2}}\,v\,\psi\Big) \;=\; & \frac{1}{R}\,v\,[\psi,u]^{st} \\
&- \frac{1}{R^{2}}\,v\,F_0\,\partial_{\phi} u \\
&+ \frac{1}{R^{2}}\,v\,\eta(\check{T}_e, \rho, \rho_{imp})\,\left(j - j_{src} - j_b\right) \\
&+ \frac{1}{R}\,\eta_{num}(\check{T}_e)\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} j\right) \\
&- \frac{1}{R^{2}}\,v\,\eta(\check{T}_e, \rho, \rho_{imp})\,j_{RE}^{ind} \\
&- \frac{2}{\rho^{\mathrm{c}}\,B^2(\psi)\,R^{3}}\,v\,\tau_{IC}\,F_0^{2}\,[\psi,p_e]^{st} \\
&+ \frac{2}{\rho^{\mathrm{c}}\,B^2(\psi)\,R^{4}}\,v\,\tau_{IC}\,F_0^{3}\,\partial_{\phi} p_e
\end{aligned}
$$

[Definitions and term groups](#details-psi)

<a id="eq-u"></a>

### Perpendicular momentum equation (`var_u`)

$$
\begin{aligned}
R^{2}\,\left(-\rho^{\mathrm{c}}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \partial_t u\right) - f_{cons}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right)\,\partial_t \rho\right) \;=\; &- \frac{1}{2}\,R\,\left(\left(\partial_{R} u\right)^{2} + \left(\partial_{Z} u\right)^{2}\right)\,[v,\hat{\rho}] \\
&- R^{3}\,\rho\,\omega\,[v,u]^{st} \\
&+ \frac{1}{R}\,v\,[\psi,j]^{st} \\
&- \frac{1}{R^{2}}\,v\,F_0\,\partial_{\phi} j \\
&+ R\,[v,p]^{st} \\
&- \mu_\perp(\check{T}_e)\,R^{2}\,f_\mu^{old}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \omega\right) \\
&- 2\,\mu_\perp(\check{T}_e)\,R\,f_\mu^{new}\,\omega\,\partial_{R} v \\
&- \mu_\perp(\check{T}_e)\,f_\mu^{new}\,\left(\partial_{R} v\,\partial_{R}\partial_{\phi}\partial_{\phi} u + \partial_{Z} v\,\partial_{Z}\partial_{\phi}\partial_{\phi} u\right) \\
&- \frac{1}{R}\,\mu_{num}(\check{T}_e)\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \omega \\
&- 2\,v\,\tau_{IC}\,R^{3}\,[p_i,\omega]^{st} \\
&- 2\,\tau_{IC}\,R^{2}\,\partial_{Z} p_i\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right) \\
&- 2\,v\,\tau_{IC}\,R^{3}\,\left(\partial_{R}\partial_{Z} u\,\left(\partial_{R}\partial_{R} p_i - \partial_{Z}\partial_{Z} p_i\right) - \partial_{R}\partial_{Z} p_i\,\left(\partial_{R}\partial_{R} u - \partial_{Z}\partial_{Z} u\right)\right) \\
&+ \mu_\perp'(\check{T}_e)\,W_{dia}\,\left(\nabla_{\mathrm{pol}} \begin{cases} T_i & \text{if } \texttt{with TiTe} \\ \frac{T}{2} & \text{otherwise} \end{cases}\cdot\nabla_{\mathrm{pol}} v\right) \\
&+ \mu_\perp(\check{T}_e)\,W_{dia}\,\nabla^2_{\mathrm{pol}} v \\
&+ f_{cons}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right)\,\left(-R\,[\hat{\rho},u] + F_0\,\left(\rho\,\partial_{\phi} v_\parallel + v_\parallel\,\partial_{\phi} \rho\right) + R\,\rho\,[v_\parallel,\psi] + R\,v_\parallel\,[\rho,\psi]\right) \\
&- v\,\left(P_\parallel^{RE} + P_\perp^{RE}\right) \\
&+ R\,\left(-\Pi_R^{RE}\,\partial_{Z} v + \Pi_Z^{RE}\,\partial_{R} v\right) \\
&+ R^{2}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right)\,\left(\left(1 - \delta_{n}\right)\,S_{neut} + \left(1 - f_{cons}\right)\,\left(S_p + S_{pellet} + S_{bg,drift} + S_{imp,drift}\right)\right) \\
&- \frac{1}{4}\,c_{TG}^{u}\,R^{4}\,\rho\,[\omega,u]\,[v,u]\,\Delta t \\
&- \frac{1}{4}\,c_{TG}^{u}\,R^{2}\,\omega\,f_{cons}\,[\hat{\rho},u]\,[v,u]\,\Delta t
\end{aligned}
$$

[Definitions and term groups](#details-u)

<a id="eq-zj"></a>

### Current definition (`var_zj`)

$$
\begin{aligned}
0 \;=\; &\frac{1}{R^{2}}\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \psi + v\,j\right)
\end{aligned}
$$

[Definitions and term groups](#details-zj)

<a id="eq-w"></a>

### Vorticity definition (`var_w`)

$$
\begin{aligned}
0 \;=\; &\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u + v\,\omega
\end{aligned}
$$

[Definitions and term groups](#details-w)

<a id="eq-rho"></a>

### Density equation (`var_rho`)

$$
\begin{aligned}
\frac{\partial}{\partial t}\Big(v\,\rho\Big) \;=\; & v\,\left(S_p + S_{pellet} + S_{bg,drift} + S_{imp,drift} + S_\rho^{aux}\right) \\
&+ v\,n_e\,\rho_n\,S_{ion}(\check{T}_e) \\
&- v\,n_e\,\rho_{main}\,S_{rec}(\check{T}_e) \\
&+ v\,\mathcal{C}_u(\rho) \\
&+ v\,\mathcal{C}_\parallel(\rho) \\
&+ \left(D_\parallel^{\mathrm{tot}} - D_\perp\right)\,\mathcal{D}_\parallel(v, \rho_{main}) \\
&- D_\perp\,\mathcal{D}_{tot}(v, \rho_{main}) \\
&+ \left(D_{\parallel,imp}^{\mathrm{tot}} - D_{\perp,imp}\right)\,\mathcal{D}_\parallel(v, \rho_{imp}) \\
&- D_{\perp,imp}\,\mathcal{D}_{tot}(v, \rho_{imp}) \\
&+ c_{TG}^{\rho}\,\mathcal{T}_n(v, \rho) \\
&- D_{\perp,num}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \rho \\
&+ 4\,v\,\tau_{IC}\,\partial_{Z} p_i^{dia} \\
&- \frac{1}{\sqrt{\left(\partial_{R} \psi\right)^{2} + \left(\partial_{Z} \psi\right)^{2}}}\,V_{pinch}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \psi\right)\,\rho
\end{aligned}
$$

[Definitions and term groups](#details-rho)

<a id="eq-vpar"></a>

### Parallel velocity equation (`var_vpar`)

$$
\begin{aligned}
v\,\left(\rho^{\mathrm{c}}\,B^2(\psi)\,\partial_t v_\parallel + \frac{1}{2}\,\rho^{\mathrm{c}}\,v_\parallel\,\partial_t B^2_{pol}(\psi) + f_{cons}\,v_\parallel\,B^2(\psi)\,\partial_t \rho\right) \;=\; &- v\,\left(\mathbf{B}\cdot\nabla p\right) \\
&+ \frac{1}{2}\,v_\parallel^{2}\,B^2\,\left(\rho\,\left(\mathbf{B}\cdot\nabla v\right) + v\,\left(\mathbf{B}\cdot\nabla \rho\right)\right) \\
&+ f_{cons}\,v\,v_\parallel\,B^2\,\left(\frac{[\hat{\rho},u]}{R} - \mathbf{B}\cdot\nabla \left(\rho\,v_\parallel\right)\right) \\
&- \mu_{\parallel,num}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} v_\parallel \\
&- \frac{1}{R^{2}\,B^2}\,\mu_{\parallel\parallel}\,F_0^{2}\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)\,\left(\mathbf{B}\cdot\nabla v\right) \\
&- \mu_\parallel^{\mathrm{eff}}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}}(v_\parallel - V_{rot})\right) \\
&+ v\,M_\parallel^{aux} \\
&- v\,S_n\,v_\parallel\,B^2\,\left(1 - f_{cons}\right) \\
&- v\,S_\rho^{aux}\,v_\parallel\,B^2\,\left(1 - f_{cons}\right) \\
&+ \left(1 - \delta_{n}\right)\,v\,n_e\,v_\parallel\,B^2\,\left(\left(\rho - \rho_{imp}\right)\,S_{rec}(\check{T}_e) - \rho_n\,S_{ion}(\check{T}_e)\right) \\
&- \frac{1}{4}\,c_{TG}^{v_\parallel}\,\rho\,v_\parallel^{2}\,B^2\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)\,\left(\mathbf{B}\cdot\nabla v\right)\,\Delta t \\
&- \frac{1}{4}\,c_{TG}^{v_\parallel}\,v\,v_\parallel^{2}\,B^2\,\left(1 - f_{cons}\right)\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)\,\left(\mathbf{B}\cdot\nabla \rho\right)\,\Delta t \\
&- \frac{1}{4}\,c_{TG}^{v_\parallel}\,v_\parallel^{3}\,B^2\,f_{cons}\,\left(\mathbf{B}\cdot\nabla \rho\right)\,\left(\mathbf{B}\cdot\nabla v\right)\,\Delta t \\
&+ \frac{1}{\sqrt{\left(\partial_{R} \psi\right)^{2} + \left(\partial_{Z} \psi\right)^{2}}}\,V_{pinch}\,\left(\nabla_{\mathrm{pol}} \psi\cdot\nabla_{\mathrm{pol}} v_\parallel\right)\,\rho\,v
\end{aligned}
$$

[Definitions and term groups](#details-vpar)

<a id="eq-rhoimp"></a>

### Impurity density equation (`var_rhoimp`)

$$
\begin{aligned}
\frac{\partial}{\partial t}\Big(v\,\rho_{imp}\Big) \;=\; & v\,S_{imp,drift} \\
&+ v\,\mathcal{C}_u(\rho_{imp}) \\
&+ v\,\mathcal{C}_\parallel(\rho_{imp}) \\
&+ \left(D_{\parallel,imp}^{\mathrm{tot}} - D_{\perp,imp}\right)\,\mathcal{D}_\parallel(v, \rho_{imp}) \\
&- D_{\perp,imp}\,\mathcal{D}_{tot}(v, \rho_{imp}) \\
&+ c_{TG}^{\rho_{imp}}\,\mathcal{T}_n(v, \rho_{imp}) \\
&- D_{\perp,num}^{n}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \rho_{imp}
\end{aligned}
$$

[Definitions and term groups](#details-rhoimp)

<a id="eq-rhon"></a>

### Neutral density equation (`var_rhon`)

$$
\begin{aligned}
\frac{\partial}{\partial t}\Big(v\,\rho_n\Big) \;=\; &- D_{n,R}\,\partial_{R} v\,\partial_{R} \rho_n \\
&- D_{n,Z}\,\partial_{Z} v\,\partial_{Z} \rho_n \\
&- \frac{1}{R^{2}}\,D_{n,\phi}\,\partial_{\phi} v\,\partial_{\phi} \rho_n \\
&+ \delta_{n}\,v\,\left(\mathcal{C}_u(\rho_n) + \mathcal{C}_\parallel(\rho_n)\right) \\
&- v\,n_e\,\rho_n\,S_{ion}(\check{T}_e) \\
&+ v\,n_e\,\rho_{main}\,S_{rec}(\check{T}_e) \\
&+ v\,S_n^{drift} \\
&- D_{\perp,num}^{n}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \rho_n
\end{aligned}
$$

[Definitions and term groups](#details-rhon)

<a id="eq-Ti"></a>

### Ion energy equation (`var_Ti`)

$$
\begin{aligned}
\frac{\partial}{\partial t}\Big(v\,\left(\rho^{\mathrm{c}} + \alpha_i\,\rho_{imp}^{\mathrm{c}}\right)\,T_i\Big) \;=\; & v\,\mathcal{C}_u^{p}(p_i) \\
&+ v\,\mathcal{C}_\parallel^{p}(p_i) \\
&+ \left(\kappa_{\parallel,i}(T_i) - \kappa_{\perp,i}(\rho)\right)\,\mathcal{D}_\parallel(v, T_i) \\
&- \kappa_{\perp,i}(\rho)\,\mathcal{D}_{tot}(v, T_i) \\
&+ v\,\left(H_i + H_i^{aux}\right) \\
&+ v\,Q_{ie}(T_i, T_e, \rho, \rho_{imp}) \\
&- v\,T_i\,\left(\rho^{\mathrm{c}}\right)^{2}\,S_{rec}(T_e) \\
&+ \frac{1}{2}\,v\,\left(\Gamma - 1\right)\,\left(v_\parallel^{2}\,B^2(\psi) + R^{2}\,\left(\left(\partial_{R} u\right)^{2} + \left(\partial_{Z} u\right)^{2}\right)\right)\,S_{kin}(T_e) \\
&+ Q_{\mu\parallel}(v) \\
&+ Q_{\mu\perp}(v, \mu_{heat}(T_e)) \\
&+ Q_{kin}(v) \\
&+ \mathcal{T}_p(v, p_i; c_{TG}^{T_i}) \\
&- \kappa_{\perp,num}^{i}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} T_i \\
&+ \mathcal{H}(v, \epsilon_i^{\mathrm{floor}}(T_i), T_i^{\mathrm{floor}}(T_i), T_{\mathrm{min}})
\end{aligned}
$$

[Definitions and term groups](#details-Ti)

<a id="eq-Te"></a>

### Electron energy equation (`var_Te`)

$$
\begin{aligned}
\frac{\partial}{\partial t}\Big(v\,\left(\rho^{\mathrm{c}}\,T_e + \rho_{imp}^{\mathrm{c}}\,A_e(T_e) + \left(\Gamma - 1\right)\,W_{ion}\right)\Big) \;=\; & v\,\mathcal{C}_u^{p}(p_e) \\
&+ v\,\mathcal{C}_\parallel^{p}(p_e) \\
&+ \left(\kappa_{\parallel,e}(T_e) - \kappa_{\perp,e}(\rho)\right)\,\mathcal{D}_\parallel(v, T_e) \\
&- \kappa_{\perp,e}(\rho)\,\mathcal{D}_{tot}(v, T_e) \\
&+ v\,\left(H_e + H_e^{aux} + P_{teleport}\right) \\
&- v\,\xi_{ion}\,\left(\rho + \alpha_e(T_e)\,\rho_{imp}\right)\,\rho_n\,S_{ion}(T_e) \\
&+ v\,Q_{ei}(T_i, T_e, \rho, \rho_{imp}) \\
&+ v\,\left(\Gamma - 1\right)\,\eta_{ohm}(T_e, \rho, \rho_{imp})\,\left(\frac{1}{R}\left(j - j_{RE}\right)\right)^{2} \\
&- v\,n_e^{\mathrm{c}}\,\rho_n^{\mathrm{c}}\,L_{rays}(T_e) \\
&- v\,n_e^{\mathrm{c}}\,\left(\rho^{\mathrm{c}} - \rho_{imp}^{\mathrm{c}}\right)\,L_{cont}(T_e) \\
&- v\,n_e^{\mathrm{c}}\,f_{rad}^{bg}(T_e) \\
&- v\,n_e^{\mathrm{c}}\,\rho_{imp}^{\mathrm{c}}\,L_{rad}(T_e) \\
&+ \mathcal{Q}_{ion}(v, W_{ion}, T_e) \\
&+ \mathcal{T}_p(v, p_e; c_{TG}^{T_e}) \\
&- \kappa_{\perp,num}^{e}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} T_e \\
&+ \mathcal{H}(v, \epsilon_e^{\mathrm{floor}}(T_e), T_e^{\mathrm{floor}}(T_e), T_{\mathrm{min}})
\end{aligned}
$$

[Definitions and term groups](#details-Te)

<a id="eq-T"></a>

### Single-temperature energy equation (`var_T`)

$$
\begin{aligned}
\frac{\partial}{\partial t}\Big(v\,\left(\rho^{\mathrm{c}}\,T + \rho_{imp}^{\mathrm{c}}\,A_{imp}(T) + \left(\Gamma - 1\right)\,W_{ion}\right)\Big) \;=\; & v\,\mathcal{C}_u^{p}(p) \\
&+ v\,\mathcal{C}_\parallel^{p}(p) \\
&+ \left(\kappa_\parallel(T) - \kappa_\perp(\rho)\right)\,\mathcal{D}_\parallel(v, T) \\
&- \kappa_\perp(\rho)\,\mathcal{D}_{tot}(v, T) \\
&+ v\,\left(H + H^{aux} + P_{teleport}\right) \\
&- v\,\xi_{ion}\,\left(\rho + \alpha_e(T)\,\rho_{imp}\right)\,\rho_n\,S_{ion}(T) \\
&- \frac{1}{2}\,\left(\Gamma - 1\right)\,v\,T\,\left(\rho^{\mathrm{c}}\right)^{2}\,S_{rec}(T) \\
&+ v\,\left(\Gamma - 1\right)\,\eta_{ohm}(T, \rho, \rho_{imp})\,\left(\frac{1}{R}\left(j - j_{RE}\right)\right)^{2} \\
&- v\,n_e^{\mathrm{c}}\,\rho_n^{\mathrm{c}}\,L_{rays}(T) \\
&- v\,n_e^{\mathrm{c}}\,\left(\rho^{\mathrm{c}} - \rho_{imp}^{\mathrm{c}}\right)\,L_{cont}(T) \\
&- v\,n_e^{\mathrm{c}}\,f_{rad}^{bg}(T) \\
&- v\,n_e^{\mathrm{c}}\,\rho_{imp}^{\mathrm{c}}\,L_{rad}(T) \\
&+ \mathcal{Q}_{ion}(v, W_{ion}, T) \\
&+ \frac{1}{2}\,v\,\left(\Gamma - 1\right)\,\left(v_\parallel^{2}\,B^2(\psi) + R^{2}\,\left(\left(\partial_{R} u\right)^{2} + \left(\partial_{Z} u\right)^{2}\right)\right)\,S_{kin}(T) \\
&+ Q_{\mu\parallel}(v) \\
&+ Q_{\mu\perp}(v, \mu_{heat}(T)) \\
&+ Q_{kin}(v) \\
&+ \mathcal{T}_p(v, p; c_{TG}^{T}) \\
&- \kappa_{\perp,num}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} T \\
&+ \mathcal{H}(v, \epsilon^{\mathrm{floor}}(T), T^{\mathrm{floor}}(T), T_{\mathrm{min}})
\end{aligned}
$$

[Definitions and term groups](#details-T)


<!-- END GENERATED: equations -->

## Status and options

### Status

| Item | Equations | Status |
|---|---|---|
| Differences with the Fortran | several | **Recorded**, see [Known differences](#known-differences). |
| Neoclassical terms | $u$, $v_\parallel$ | **Not implemented.** The Fortran neoclassical terms (profiles `amu_neo_prof`, `aki_neo_prof`) are not transcribed, and the checker does not audit the NEO branch. The `include_neo` argument of `momentum_equation_2` has no effect. |
| Inward pinch | $\rho$, $v_\parallel$ | **Under review.** Whether the $v_\parallel$ pinch term belongs there is open, and the Fortran Jacobian of both pinch terms disagrees with the checker (findings 4 and 9). |
| Momentum of the neutral sources | $v_\parallel$ | **Suspected defect, to confirm.** Probably needs the factor $(1-f_{cons})$ that the other particle-source momentum terms carry. |
| Freezing $u$ in the time derivative | $u$ | **To confirm.** The element routine holds $u$ at its current value in the time derivative of the conservative-form term. |
| Impurity pressure in the $\tau_{IC}$ tangents | $\psi$, $u$, $\rho$ | **Accepted.** The residual is complete; the Fortran Jacobian omits the impurity part (findings 1 and 2). Impurities and $\tau_{IC} \neq 0$ are not used together. |

#### Known differences

The findings still recorded in the checker's reference
(`util/equation_codegen/reference/model600_discrepancies.json`), and the
blocks of the element routine they affect. A finding with residual blocks
changes the converged solution; one with only Jacobian blocks affects only
the Newton convergence. The table is generated from the reference.

<!-- BEGIN GENERATED: findings -->

| Finding | Residual blocks | Jacobian blocks |
|---:|---|---|
| [1](https://github.com/iterorganization/JOREK/blob/develop/util/equation_codegen/JOREK_FINDINGS.md) | none | `amat(var_rho,var_rhoimp)`, `amat(var_rho,var_t)`, `amat(var_rho,var_ti)`, `amat(var_u,var_rhoimp)`, `amat(var_u,var_t)`, `amat(var_u,var_ti)` |
| [2](https://github.com/iterorganization/JOREK/blob/develop/util/equation_codegen/JOREK_FINDINGS.md) | none | `amat(var_psi,var_rhoimp)`, `amat(var_psi,var_t)`, `amat(var_psi,var_te)`, `amat_n(var_psi,var_rhoimp)`, `amat_n(var_psi,var_t)`, `amat_n(var_psi,var_te)` |
| [4](https://github.com/iterorganization/JOREK/blob/develop/util/equation_codegen/JOREK_FINDINGS.md) | none | `amat(var_rho,var_psi)`, `amat(var_vpar,var_psi)` |
| [9](https://github.com/iterorganization/JOREK/blob/develop/util/equation_codegen/JOREK_FINDINGS.md) | none | `amat(var_vpar,var_rho)`, `amat(var_vpar,var_vpar)` |

<!-- END GENERATED: findings -->

### Options

The compile-time switches of `mod_model_settings.f90` (see
[Variables](#variables)) are genuine branches: they remove a variable and its
equation, and in the element routine they gate whole blocks of terms. All other
features are switched by the coefficients below, as they are named in the
element routine (the input parameters that set them are documented with the
model); a disabled feature has a zero coefficient, and its terms are still
assembled.

| Feature | Controlled by | Mechanism |
|---|---|---|
| Two-temperature model | `with_TiTe` | branch: equations for $T_i$, $T_e$ instead of $T$ |
| Density, parallel velocity, neutrals, impurities | `with_rho`, `with_vpar`, `with_neutrals`, `with_impurities` | branch: the variable is replaced by a constant and its equation is not assembled |
| Diamagnetic terms | $\tau_{IC}$ (`tauIC`) | coefficient |
| Conservative form of the momentum equations | $f_{cons}$ (`fact_conservative_u`) | coefficient, $0$ or $1$ |
| Convection of the neutrals with the plasma flow, and the momentum they carry | $\delta_n$ (`delta_n_convection`) | the neutral equation convects $\rho_n$ with $\delta_n$ (on for $\delta_n = 1$); the momentum of the neutral sources carries $(1-\delta_n)$ (on for $\delta_n = 0$) |
| Inward pinch | $V_{pinch}$ (`V_prof_pinch`) | profile, zero when off |
| Taylor-Galerkin stabilization | $c_{TG}$ (`tgnum_*`) | coefficient per equation |
| Numerical hyper-diffusion | $D_{\perp,num}$, $\kappa_{\perp,num}$, $\eta_{num}$, $\mu_{num}$, $\mu_{\parallel,num}$ | coefficients |
| Shock capturing | $D_{\parallel,sc}\,\tau_{sc}$, $D_{\parallel,imp,sc}\,\tau_{sc}$, $\mu_{\parallel,sc}\,\tau_{sc}$ | products of a coefficient and the shock indicator $\tau_{sc}$ |
| Runaway-electron coupling | $j_{RE}$, $P^{RE}$, $\Pi^{RE}$ (`aux_jre`, `aux_P_*_re`, `aux_divPI*_perp`) | coupling values, zero without runaways |
| Kinetic neutral and impurity coupling | $S_\rho^{aux}$, $M_\parallel^{aux}$, $H^{aux}$ (`aux_rho0`, `aux_mom_par0`, `aux_E0*`) | coupling values, zero without the kinetic model |
| Neoclassical terms | `amu_neo_prof`, `aki_neo_prof` in the Fortran | not transcribed (see Status) |

## Conventions

The weak form assumes:

- $v$ is a finite-element basis function in space (of the element
  coordinates and of $\phi$); it does not depend on time.
- The fields are periodic in $\phi$.
- The grid, hence $R$ and $\mathcal{J}$, does not change in time.
- Integration by parts moves derivatives onto $v$; the boundary terms it
  produces are not part of these volume terms. Boundary conditions are imposed
  separately (`models/model600/mod_boundary_conditions.f90`).

**Time derivatives.** The equations at the top of the page write every time
derivative out. In the momentum and parallel velocity equations it is not the
derivative of a single quantity: the parallel velocity equation, for example,
is the projection of the momentum equation on $\mathbf{B}$ and has
$\rho^{\mathrm{c}}\big(B^2\,\partial_t v_\parallel + \frac{1}{2}v_\parallel\,\partial_t B^2_{pol}\big)$.
In `model600.py` such a term is written as a mass functional $A$ with
frozen factors (`freeze(x)`, shown as $\overline{x}$), and $\partial_t$ acts only on the factors without an overline:
$\partial_t(\overline{a}\,b)$ means $a\,\partial_t b$. The symbolic model never
uses the value of $A$ itself, only this derivative, and the checker derives the
Jacobian from it with the same rule. Where the density correction is inactive
($\rho^{\mathrm{c}} = \rho$), the momentum term
$-\rho\,\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}}\partial_t u - f_{cons}\,\partial_t\rho\,\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u$
is the time derivative of $\rho\,\nabla_{\mathrm{pol}} u$ for $f_{cons} = 1$
(conservative form), and $\rho\,\partial_t\nabla_{\mathrm{pol}} u$ for
$f_{cons} = 0$. Overlines appear only in these two mass functionals.

How the other terms are linearized for the Newton iteration follows from the
equations, apart from the state-dependent functions: the table
[Linearization of the state-dependent functions](#linearization-of-the-state-dependent-functions)
says which derivatives the Jacobian uses for each of them.

In `model600.py` an evolution equation is written as a pair $A$, $B$ meaning
$\frac{\partial}{\partial t}\int A = \int B$, and a constraint as $C$ meaning
$\int C = 0$. There $A$, $B$ and $C$ include the volume weight
$w_V = R\,\mathcal{J}$ (the Python name is `dV`), where $\mathcal{J}$ is the
Jacobian of the element coordinates $(s,t_{\mathrm{el}})$; the equations at the top of the
page show them divided by $w_V$. The time discretization and the Jacobian (the
`amat` blocks) are not written here: the checker derives them by
differentiating $A$ and $B$ with respect to each variable.

- $[a,b] = \partial_R a\,\partial_Z b - \partial_Z a\,\partial_R b$ is the
  poloidal bracket.
- $[a,b]^{st} = \big(\partial_s a\,\partial_{t_{\mathrm{el}}} b - \partial_{t_{\mathrm{el}}} a\,\partial_s b\big)/\mathcal{J}$
  is the same bracket computed in the element coordinates:
  $[a,b]^{st} = [a,b]^{RZ} = [a,b]$. The element routine uses both spellings
  (in `model600.py`, `poiss_bracket_st(a, b) / xjac` and `poiss_bracket(a, b)`);
  they are equal, but the checker keeps them apart because they expand into
  different monomials. Where the element-coordinate bracket appears without
  its $1/\mathcal{J}$, it is written $\mathcal{J}\,[a,b]^{st}$.
- $\nabla_{\mathrm{pol}} = (\partial_R, \partial_Z)$ is the **poloidal** gradient: it has
  no toroidal component, and
  $\nabla_{\mathrm{pol}} a\cdot\nabla_{\mathrm{pol}} b = \partial_R a\,\partial_R b + \partial_Z a\,\partial_Z b$.
  Toroidal derivatives are
  always written explicitly as $\partial_\phi$. Likewise
  $\nabla^2_{\mathrm{pol}} f = \partial_R\partial_R f + \partial_Z\partial_Z f + \partial_R f / R$
  is the poloidal part of the Laplacian.
- $\mathbf{B}\cdot\nabla$ is the full three-dimensional derivative along the field,
  including its toroidal part (see [Parallel gradient](#parallel-gradient)).
- $\overline{x}$ in a mass functional $A$ is a factor that $\partial_t$ does
  not act on (see above).
- The two diffusion operators carry opposite signs, as in `model600.py`:
  $\mathcal{D}_\parallel(v,n) = -(\mathbf{B}\cdot\nabla v)(\mathbf{B}\cdot\nabla n)/B^2$
  is the integrated-by-parts form of $+\nabla\cdot(\mathbf{b}\,\mathbf{b}\cdot\nabla n)$,
  while $\mathcal{D}_{tot}(v,f) = \nabla v\cdot\nabla f + \ldots$ is that of $-\nabla^2 f$.
  A parallel diffusion therefore appears as $+(D_\parallel - D_\perp)\,\mathcal{D}_\parallel(v, n)$
  and a perpendicular one as $-D_\perp\,\mathcal{D}_{tot}(v, n)$.
- $\Gamma$ is the ratio of specific heats. The element routine spells it both
  `GAMMA` and `gamma`; Fortran does not distinguish them.
- `with_TiTe` selects the two-temperature model (fields $T_i$, $T_e$) instead
  of the single-temperature one (field $T$); it is both a switch of the
  element routine and an argument of the equation functions of `model600.py`. `st_form` only selects which
  spelling of a bracket the element routine uses; it does not change the
  physics. The checker needs it because the two spellings expand into
  different monomials.

## Notation

Names that are not listed are rendered from their Python spelling:
`D_par_local` becomes $D_{\mathrm{par},\mathrm{local}}$.

<!-- BEGIN GENERATED: notation -->

### Work values

| Python | Symbol |
|---|---|
| `T_or_Te` | $\check{T}_e$ |
| `Te_gen` | $\hat{T}_e$ |
| `ion_temperature` | $\tilde{T}_i$ |
| `R` | $R$ |
| `xjac` | $\mathcal{J}$ |
| `dV` | $w_V$ |
| `F0` | $F_0$ |
| `tauIC` | $\tau_{IC}$ |
| `tstep` | $\Delta t$ |
| `GAMMA` | $\Gamma$ |
| `gamma` | $\Gamma$ |
| `ne` | $n_e$ |
| `rho_main` | $\rho_{main}$ |
| `rho_corr` | $\rho^{\mathrm{c}}$ |
| `rho_hat` | $\hat{\rho}$ |
| `bb2` | $B^2$ |
| `Pe` | $p_e$ |
| `pressure` | $p$ |
| `pi` | $p_i$ |
| `W_dia` | $W_{dia}$ |
| `ion_density` | $n_i$ |
| `ion_pressure` | $p_i$ |
| `electron_pressure` | $p_e$ |
| `electron_density` | $n_e^{\mathrm{c}}$ |
| `ionization_energy` | $W_{ion}$ |
| `visco_par_eff` | $\mu_\parallel^{\mathrm{eff}}$ |
| `source_dens_tot` | $S_n$ |
| `D_par_tot` | $D_\parallel^{\mathrm{tot}}$ |
| `D_par_imp_tot` | $D_{\parallel,imp}^{\mathrm{tot}}$ |
| `d_par_tot` | $\hat{D}_\parallel$ |
| `d_par_imp_tot` | $\hat{D}_{\parallel,imp}$ |
| `alpha_e_value` | $\alpha_e$ |
| `sion_rate` | $S_{ion}$ |
| `srec_rate` | $S_{rec}$ |
| `grad_v` | $\nabla_{\mathrm{pol}} v$ |
| `grad_u` | $\nabla_{\mathrm{pol}} u$ |
| `grad_omega` | $\nabla_{\mathrm{pol}} \omega$ |
| `grad_vpar` | $\nabla_{\mathrm{pol}} v_\parallel$ |
| `lap_v` | $\nabla^2_{\mathrm{pol}} v$ |
| `lap_omega` | $\nabla^2_{\mathrm{pol}} \omega$ |
| `neutral_sources` | $S_{neut}$ |
| `current_source` | $j_{src}$ |
| `Jb` | $j_b$ |
| `visco_par` | $\mu_\parallel$ |
| `visco_par_num` | $\mu_{\parallel,num}$ |
| `visco_par_par` | $\mu_{\parallel\parallel}$ |
| `visco_par_sc_num` | $\mu_{\parallel,sc}$ |
| `D_prof` | $D_\perp$ |
| `D_prof_imp` | $D_{\perp,imp}$ |
| `D_par_local` | $D_\parallel$ |
| `D_par_local_imp` | $D_{\parallel,imp}$ |
| `D_par_sc_num` | $D_{\parallel,sc}$ |
| `D_par_imp_sc_num` | $D_{\parallel,imp,sc}$ |
| `tau_sc` | $\tau_{sc}$ |
| `V_prof_pinch` | $V_{pinch}$ |
| `fact_conservative_u` | $f_{cons}$ |
| `delta_n_convection` | $\delta_{n}$ |
| `dV_dpsi_source` | $\frac{\mathrm{d}V_{rot}}{\mathrm{d}\psi}$ |
| `particle_source` | $S_p$ |
| `source_pellet` | $S_{pellet}$ |
| `source_bg_drift` | $S_{bg,drift}$ |
| `source_imp_drift` | $S_{imp,drift}$ |
| `aux_rho0` | $S_\rho^{aux}$ |
| `aux_mom_par0` | $M_\parallel^{aux}$ |
| `heat_source_total` | $H$ |
| `heat_source_i` | $H_i$ |
| `heat_source_e` | $H_e$ |
| `aux_E0` | $H^{aux}$ |
| `aux_E0_Ti` | $H_i^{aux}$ |
| `aux_E0_Te` | $H_e^{aux}$ |
| `power_dens_teleport_ju` | $P_{teleport}$ |
| `ksi_ion_norm` | $\xi_{ion}$ |
| `implicit_heat_source` | $c_{\mathrm{floor}}$ |
| `T_min_neg` | $T_{\mathrm{min}}$ |
| `Tie_min_neg` | $T_{\mathrm{min}}$ |
| `aux_jre` | $j_{RE}$ |
| `aux_jre_ind` | $j_{RE}^{ind}$ |
| `aux_P_par_re` | $P_\parallel^{RE}$ |
| `aux_P_perp_re` | $P_\perp^{RE}$ |
| `aux_divPIR_perp` | $\Pi_R^{RE}$ |
| `aux_divPIZ_perp` | $\Pi_Z^{RE}$ |
| `visco_fact_old` | $f_\mu^{old}$ |
| `visco_fact_new` | $f_\mu^{new}$ |
| `visco_par_heating` | $\mu_{\parallel,heat}$ |
| `D_perp_num_psin` | $D_{\perp,num}$ |
| `Dn0x` | $D_{n,R}$ |
| `Dn0y` | $D_{n,Z}$ |
| `Dn0p` | $D_{n,\phi}$ |
| `source_neutral_drift` | $S_n^{drift}$ |
| `Dn_perp_num` | $D_{\perp,num}^{n}$ |
| `ZK_perp_num_psin` | $\kappa_{\perp,num}$ |
| `ZK_i_perp_num_psin` | $\kappa_{\perp,num}^{i}$ |
| `ZK_e_perp_num_psin` | $\kappa_{\perp,num}^{e}$ |
| `tgnum_u` | $c_{TG}^{u}$ |
| `tgnum_rho` | $c_{TG}^{\rho}$ |
| `tgnum_vpar` | $c_{TG}^{v_\parallel}$ |
| `tgnum_rhoimp` | $c_{TG}^{\rho_{imp}}$ |
| `tgnum_T` | $c_{TG}^{T}$ |
| `tgnum_Ti` | $c_{TG}^{T_i}$ |
| `tgnum_Te` | $c_{TG}^{T_e}$ |
| `psi_gradient` | $\nabla_{\mathrm{pol}} \psi$ |
| `rotation_shear` | $\nabla_{\mathrm{pol}}(v_\parallel - V_{rot})$ |

### State-dependent functions

A function of the state is shown with its arguments, e.g. $S_{ion}(T)$;
$x^{\mathrm{c}}$ is a density corrected for negative values.

| Python | Symbol |
|---|---|
| `eta` | $\eta$ |
| `eta_num_T` | $\eta_{num}$ |
| `visco` | $\mu_\perp$ |
| `visco_num` | $\mu_{num}$ |
| `alpha_i_state` | $\alpha_i$ |
| `dvisco_state` | $\mu_\perp'$ |
| `Sion_rate` | $S_{ion}$ |
| `Srec_rate` | $S_{rec}$ |
| `alpha_e_state` | $\alpha_e$ |
| `alpha_e_bis_state` | $\alpha_e'$ |
| `alpha_e_temperature` | $A_e$ |
| `alpha_imp_bis_state` | $\alpha_{imp}'$ |
| `alpha_imp_temperature` | $A_{imp}$ |
| `corr_neg_dens` | ${0}^{\mathrm{c}}$ |
| `corr_neg_dens_imp` | ${0}^{\mathrm{c}}$ |
| `corr_neg_dens_n` | ${0}^{\mathrm{c}}$ |
| `W_dia_single` | $W_{dia}$ |
| `W_dia_two` | $W_{dia}$ |
| `ZKi_par` | $\kappa_{\parallel,i}$ |
| `ZKi_perp` | $\kappa_{\perp,i}$ |
| `visco_heating` | $\mu_{heat}$ |
| `Ti_e_exchange` | $Q_{ie}$ |
| `Ti_floor` | $T_i^{\mathrm{floor}}$ |
| `Ti_floor_exp` | $\epsilon_i^{\mathrm{floor}}$ |
| `E_ion_bg_state` | $E_{ion}^{bg}$ |
| `E_ion_state` | $E_{ion}$ |
| `ZKe_par` | $\kappa_{\parallel,e}$ |
| `ZKe_perp` | $\kappa_{\perp,e}$ |
| `eta_ohm_e` | $\eta_{ohm}$ |
| `LradDrays` | $L_{rays}$ |
| `LradDcont` | $L_{cont}$ |
| `frad_bg_state` | $f_{rad}^{bg}$ |
| `Lrad_state` | $L_{rad}$ |
| `Te_i_exchange` | $Q_{ei}$ |
| `Te_floor` | $T_e^{\mathrm{floor}}$ |
| `Te_floor_exp` | $\epsilon_e^{\mathrm{floor}}$ |
| `ZK_par` | $\kappa_\parallel$ |
| `ZK_perp` | $\kappa_\perp$ |
| `T_floor` | $T^{\mathrm{floor}}$ |
| `T_floor_exp` | $\epsilon^{\mathrm{floor}}$ |

<!-- END GENERATED: notation -->

### Linearization of the state-dependent functions

Generated from the `external_function` declarations of `model600.py`. The
Jacobian differentiates a function through the derivatives the element
routine supplies for it; "current branch" means the derivative of the branch
(or clipping regime) selected at the current state, so a clipped constant has
zero derivative.

<!-- BEGIN GENERATED: linearization -->

| DSL | Arguments | Fortran value | Derivatives supplied | Jacobian |
|---|---|---|---|---|
| `eta` | `T`, `rho`, `rhoimp` | `eta_T` | `deta_dT` for `T`, `deta_dr0` for `rho`, `deta_drimp0` for `rhoimp` | differentiated (current branch) |
| `eta_num_T` | `T` | `eta_num_T` | `deta_num_dT` for `T` | differentiated (current branch) |
| `visco` | `T` | `visco_T` | `dvisco_dT` for `T` | differentiated (current branch) |
| `visco_num` | `T` | `visco_num_T` | none | not differentiated |
| `alpha_i_state` | none | `alpha_i` | none | not differentiated |
| `dvisco_state` | `T` | `dvisco_dT` | `d2visco_dT2` for `T` | differentiated (current branch) |
| `Sion_rate` | `T` | `Sion_T` | `dSion_dT` for `T` | differentiated (current branch) |
| `Srec_rate` | `T` | `Srec_T` | `dSrec_dT` for `T` | differentiated (current branch) |
| `alpha_e_state` | `T` | `alpha_e` | `dalpha_e_dT` for `T` | differentiated (current branch) |
| `alpha_e_bis_state` | `T` | `alpha_e_bis` | `alpha_e_tri` for `T` | differentiated (current branch) |
| `alpha_e_temperature` | `T` | `alpha_e_T` | `alpha_e_bis` for `T` | differentiated (current branch) |
| `alpha_imp_bis_state` | `T` | `alpha_imp_bis` | `alpha_imp_tri` for `T` | differentiated (current branch) |
| `alpha_imp_temperature` | `T` | `alpha_imp_T` | `alpha_imp_bis` for `T` | differentiated (current branch) |
| `corr_neg_dens` | `rho` | `corr_neg_dens` | `dr0_corr_dn` for `rho` | differentiated (current branch) |
| `W_dia_single` | `rho`, `T` | `W_dia` | `W_dia_rho` for `rho`, `W_dia_T` for `T` | differentiated (current branch) |
| `W_dia_two` | `rho`, `Ti` | `W_dia` | `W_dia_rho` for `rho`, `W_dia_Ti` for `Ti` | differentiated (current branch) |
| `ZKi_par` | `Ti` | `ZKi_par_T` | `dZKi_par_dT` for `Ti` | differentiated (current branch) |
| `ZKi_perp` | `rho` | `ZKi_prof` | `dZKi_prof_drho` for `rho` | differentiated (current branch) |
| `visco_heating` | `Te` | `visco_T_heating` | `dvisco_dT_heating` for `Te` | differentiated (current branch) |
| `Ti_e_exchange` | `Ti`, `Te`, `rho`, `rhoimp` | `dTi_e` | `ddTi_e_dTi` for `Ti`, `ddTi_e_dTe` for `Te`, `ddTi_e_drho` for `rho`, `ddTi_e_drhoimp` for `rhoimp` | differentiated (current branch) |
| `Ti_floor` | `Ti` | `Ti0_floor` | `dTi_floor` for `Ti` | differentiated (current branch) |
| `Ti_floor_exp` | `Ti` | `Ti_floor_exp` | `dTi_floor_exp` for `Ti` | differentiated (current branch) |
| `corr_neg_dens_imp` | `rhoimp` | `corr_neg_dens_imp` | `drimp0_corr_dn` for `rhoimp` | differentiated (current branch) |
| `E_ion_bg_state` | none | `E_ion_bg` | none | not differentiated |
| `ZKe_par` | `Te` | `ZKe_par_T` | `dZKe_par_dT` for `Te` | differentiated (current branch) |
| `ZKe_perp` | `rho` | `ZKe_prof` | `dZKe_prof_drho` for `rho` | differentiated (current branch) |
| `eta_ohm_e` | `Te`, `rho`, `rhoimp` | `eta_T_ohm` | `deta_dT_ohm` for `Te`, `deta_dr0_ohm` for `rho`, `deta_drimp0_ohm` for `rhoimp` | differentiated (current branch) |
| `dE_ion_dT_state` | `Te` | `dE_ion_dT` | none | not differentiated |
| `E_ion_state` | `Te` | `E_ion` | `dE_ion_dT` for `Te` | differentiated (current branch) |
| `LradDrays` | `Te` | `LradDrays_T` | `dLradDrays_dT` for `Te` | differentiated (current branch) |
| `LradDcont` | `Te` | `LradDcont_corr` | `dLradDcont_dT_corr` for `Te` | differentiated (current branch) |
| `frad_bg_state` | `Te` | `frad_bg` | `dfrad_bg_dT` for `Te` | differentiated (current branch) |
| `Lrad_state` | `Te` | `Lrad` | `dLrad_dT` for `Te` | differentiated (current branch) |
| `Te_i_exchange` | `Ti`, `Te`, `rho`, `rhoimp` | `dTe_i` | `ddTe_i_dTi` for `Ti`, `ddTe_i_dTe` for `Te`, `ddTe_i_drho` for `rho`, `ddTe_i_drhoimp` for `rhoimp` | differentiated (current branch) |
| `Te_floor` | `Te` | `Te0_floor` | `dTe_floor` for `Te` | differentiated (current branch) |
| `Te_floor_exp` | `Te` | `Te_floor_exp` | `dTe_floor_exp` for `Te` | differentiated (current branch) |
| `corr_neg_dens_n` | `rhon` | `corr_neg_dens_n` | `drn0_corr_dn` for `rhon` | differentiated (current branch) |
| `ZK_par` | `T` | `ZK_par_T` | `dZK_par_dT` for `T` | differentiated (current branch) |
| `ZK_perp` | `rho` | `ZK_prof` | `dZK_prof_drho` for `rho` | differentiated (current branch) |
| `T_floor` | `T` | `T0_floor` | `dT_floor` for `T` | differentiated (current branch) |
| `T_floor_exp` | `T` | `T_floor_exp` | `dT_floor_exp` for `T` | differentiated (current branch) |

<!-- END GENERATED: linearization -->

## Shared operators

Term groups that appear in several equations, grouped by topic. The symbol on
the left of each definition is how the group is written in the equations
above; $f$, $n$ and $p$ stand for any scalar, density and pressure. Each is a
`_helper` function of `model600.py`.

<!-- BEGIN GENERATED: operators -->

### Magnetic field

#### Magnetic field strength

The square of the total field, with $\mathbf{B} = F_0\,\nabla\phi + \nabla\psi\times\nabla\phi$.

$$
\begin{aligned}
B^2(\psi) &= \frac{1}{R^{2}}\left(F_0^{2} + \left(\partial_{R} \psi\right)^{2} + \left(\partial_{Z} \psi\right)^{2}\right)
\end{aligned}
$$

Source: `_B2`.

#### Poloidal magnetic field strength

The square of the poloidal field, $\lvert\nabla_{\mathrm{pol}}\psi\times\nabla\phi\rvert^2$. Since $F_0$ is constant, $\partial_t B^2 = \partial_t B^2_{pol}$.

$$
\begin{aligned}
B^2_{pol}(\psi) &= \frac{1}{R^{2}}\left(\left(\partial_{R} \psi\right)^{2} + \left(\partial_{Z} \psi\right)^{2}\right)
\end{aligned}
$$

Source: `_B2_pol`.

#### Parallel gradient

The derivative along the field. With `st_form` the bracket is written in element coordinates, $[f,\psi]^{st}$, which is the same quantity.

$$
\begin{aligned}
\mathbf{B}\cdot\nabla f &= \frac{1}{R}\left(\frac{1}{R}\,F_0\,\partial_{\phi} f + [f,\psi]\right)
\end{aligned}
$$

Source: `_B_dot_grad`.

#### Weighted density

The density weighted by $R^2$, as it enters the momentum equation.

$$
\begin{aligned}
\hat{\rho} &= R^{2}\,\rho
\end{aligned}
$$

Source: `_rho_hat`.

### Convection

#### Convection by the ExB flow

Advection and compression of a density by the $E\times B$ flow of stream function $u$.

$$
\begin{aligned}
\mathcal{C}_u(n) &= R\,[n,u] + 2\,n\,\partial_{Z} u
\end{aligned}
$$

Source: `_u_convection`.

#### Parallel convection

Advection and compression of a density by the parallel flow $v_\parallel$.

$$
\begin{aligned}
\mathcal{C}_\parallel(n) &= -v_\parallel\,\left(\mathbf{B}\cdot\nabla n\right) - n\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)
\end{aligned}
$$

Source: `_parallel_convection`.

#### Pressure convection by the ExB flow

The same for a pressure: the compression carries the factor $\Gamma$.

$$
\begin{aligned}
\mathcal{C}_u^{p}(p) &= R\,[p,u] + 2\,\Gamma\,p\,\partial_{Z} u
\end{aligned}
$$

Source: `_press_u_convection`.

#### Parallel pressure convection

The same for a pressure along the field.

$$
\begin{aligned}
\mathcal{C}_\parallel^{p}(p) &= -v_\parallel\,\left(\mathbf{B}\cdot\nabla p\right) - \Gamma\,p\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)
\end{aligned}
$$

Source: `_press_parallel_convection`.

### Diffusion

#### Parallel diffusion

Diffusion along the magnetic field, after integration by parts.

$$
\begin{aligned}
\mathcal{D}_\parallel(v, n) &= -\frac{1}{B^2(\psi)}\,\left(\mathbf{B}\cdot\nabla v\right)\,\left(\mathbf{B}\cdot\nabla n\right)
\end{aligned}
$$

Source: `_par_diff_intg_by_parts`.

#### Total diffusion

Isotropic diffusion, including the toroidal derivative, after integration by parts.

$$
\begin{aligned}
\mathcal{D}_{tot}(v, f) &= \nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} f + \frac{1}{R^{2}}\,\partial_{\phi} v\,\partial_{\phi} f
\end{aligned}
$$

Source: `_diffusion_tot_intg_by_parts`.

### Numerical stabilization

#### Taylor-Galerkin stabilization of a density

Taylor-Galerkin stabilization of the convection by the $E\times B$ and parallel flows.

$$
\begin{aligned}
\mathcal{T}_n(v, n) &= -\frac{1}{4}\,R^{2}\,[n,u]\,[v,u]\,\Delta t - \frac{1}{4}\,v_\parallel^{2}\,\left(\mathbf{B}\cdot\nabla n\right)\,\left(\mathbf{B}\cdot\nabla v\right)\,\Delta t
\end{aligned}
$$

Source: `_dens_tgnum_intg_by_parts`.

#### Taylor-Galerkin stabilization of an energy equation

The same for a pressure, with the coefficient $c$ of the equation.

$$
\begin{aligned}
\mathcal{T}_p(v, p; c) &= -\frac{1}{4}\,c\,R^{2}\,[p,u]\,[v,u]\,\Delta t - \frac{1}{4}\,c\,v_\parallel^{2}\,\left(\mathbf{B}\cdot\nabla p\right)\,\left(\mathbf{B}\cdot\nabla v\right)\,\Delta t
\end{aligned}
$$

Source: `_energy_tgnum`.

#### Heating at the temperature floor

Implicit heat source that keeps a temperature away from its floor $T_{\mathrm{min}}$. $\epsilon$ and $T^{\mathrm{floor}}$ are the floor functions of the element routine (`Ti_floor_exp`, `Ti0_floor` and their electron and single-temperature counterparts).

$$
\begin{aligned}
\mathcal{H}(v, \epsilon, T^{\mathrm{floor}}, T_{\mathrm{min}}) &= c_{\mathrm{floor}}\,\left(\Gamma - 1\right)\,v\,\left(\frac{1}{2}\,T_{\mathrm{min}}\,\left(1 + \epsilon\right) - T^{\mathrm{floor}}\right)
\end{aligned}
$$

Source: `_heating_floor`.

### Diamagnetic terms

#### Diamagnetic ion pressure

The ion pressure of the diamagnetic terms, as built in `construct_pressure`. The single-temperature model evolves the total temperature, so the ion temperature is $T/2$ there. $\alpha_i$ does not depend on the temperature.

$$
\begin{aligned}
\tilde{T}_i &= \begin{cases} T_i & \text{if } \texttt{with TiTe} \\ \frac{T}{2} & \text{otherwise} \end{cases} \\
p_i^{dia} &= \left(\rho + \rho_{imp}\,\alpha_i\right)\,\tilde{T}_i
\end{aligned}
$$

Source: `_diamagnetic_pressure`.

#### Diamagnetic viscosity

The element routine's `W_dia`. It is built from the ion pressure alone, so it has no explicit electron-temperature dependence.

$$
\begin{aligned}
W_{dia} &= \begin{cases} W_{dia}(\rho, T_i) & \text{if } \texttt{with TiTe} \\ W_{dia}(\rho, T) & \text{otherwise} \end{cases}
\end{aligned}
$$

Source: `_diamagnetic_viscosity`.

### Sources and heating

#### Particle sources releasing kinetic energy

Particle sources whose kinetic energy is released into the ions: ionization of neutrals and the external sources.

$$
\begin{aligned}
S_{kin}(T) &= \left(\rho + \alpha_e(T)\,\rho_{imp}\right)\,\rho_n\,S_{ion}(T) + S_p + S_{pellet} + S_{bg,drift} + S_{imp,drift}
\end{aligned}
$$

Source: `_released_kinetic_energy`.

#### Parallel viscous heating

Heating by the parallel viscosity.

$$
\begin{aligned}
Q_{\mu\parallel}(v) &= \left(\Gamma - 1\right)\,\mu_{\parallel,heat}\,\left(v\,\left(\nabla_{\mathrm{pol}} v_\parallel\cdot\nabla_{\mathrm{pol}} v_\parallel\right) + v_\parallel\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} v_\parallel\right)\right)
\end{aligned}
$$

Source: `_parallel_viscous_heating`.

#### Perpendicular viscous heating

Heating by the perpendicular viscosity of the $u$ flow, with the heating viscosity $\mu_{heat}$.

$$
\begin{aligned}
Q_{\mu\perp}(v, \mu_{heat}) &= -\left(\Gamma - 1\right)\,v\,\mu_{heat}\,R^{2}\,f_\mu^{old}\,\left(\nabla_{\mathrm{pol}} u\cdot\nabla_{\mathrm{pol}} \omega\right) - 2\,\left(\Gamma - 1\right)\,v\,\mu_{heat}\,R\,f_\mu^{new}\,\omega\,\partial_{R} u - \left(\Gamma - 1\right)\,v\,\mu_{heat}\,f_\mu^{new}\,\left(\partial_{R} u\,\partial_{R}\partial_{\phi}\partial_{\phi} u + \partial_{Z} u\,\partial_{Z}\partial_{\phi}\partial_{\phi} u\right)
\end{aligned}
$$

Source: `_u_viscous_heating`.

#### Coupling to the kinetic neutral and impurity model

Energy and momentum handed over by the kinetic neutral and impurity model.

$$
\begin{aligned}
Q_{kin}(v) &= \frac{1}{2}\,\left(\Gamma - 1\right)\,v\,S_\rho^{aux}\,v_\parallel^{2}\,B^2(\psi) - \left(\Gamma - 1\right)\,v\,M_\parallel^{aux}\,v_\parallel
\end{aligned}
$$

Source: `_kinetic_coupling`.

#### Transport of the ionization potential energy

The ionization potential energy $W$ is convected like a density, and diffuses with the impurity and main-ion densities that carry it.

$$
\begin{aligned}
\hat{D}_\parallel &= D_\parallel + D_{\parallel,sc}\,\tau_{sc} - D_\perp \\
\hat{D}_{\parallel,imp} &= D_{\parallel,imp} + D_{\parallel,imp,sc}\,\tau_{sc} - D_{\perp,imp} \\
B^2 &= B^2(\psi) \\
\mathcal{Q}_{ion}(v, W, T) &= \left(\Gamma - 1\right)\Big(v\,R\,[W,u]^{st} \\
&\qquad + 2\,v\,W\,\partial_{Z} u \\
&\qquad - \frac{1}{R^{2}}\,v\,F_0\,v_\parallel\,\partial_{\phi} W \\
&\qquad - \frac{1}{R}\,v\,v_\parallel\,[W,\psi]^{st} \\
&\qquad - v\,W\,\left(\mathbf{B}\cdot\nabla v_\parallel\right) \\
&\qquad - \frac{1}{B^2}\,E_{ion}(T)\,\hat{D}_{\parallel,imp}\,\left(\mathbf{B}\cdot\nabla v\right)\,\left(\mathbf{B}\cdot\nabla \rho_{imp}\right) \\
&\qquad - E_{ion}(T)\,D_{\perp,imp}\,\mathcal{D}_{tot}(v, \rho_{imp}) \\
&\qquad - \frac{1}{B^2}\,E_{ion}^{bg}\,\hat{D}_\parallel\,\left(\mathbf{B}\cdot\nabla v\right)\,\left(\mathbf{B}\cdot\nabla \left(\rho - \rho_{imp}\right)\right) \\
&\qquad - E_{ion}^{bg}\,D_\perp\,\mathcal{D}_{tot}(v, \rho - \rho_{imp})\Big)
\end{aligned}
$$

Source: `_ionization_energy_transport`.


<!-- END GENERATED: operators -->

## Equation details

The same equations, group by group, as they are written in `model600.py`: the
local definitions, the mass functional $A$ and the right-hand side $B$ (or the
constraint $C$), which include the volume weight $w_V$. Each term group is a
commented block of the corresponding Python function. Inside a mass
functional, $\overline{x}$ marks a factor that $\partial_t$ does not act on
(see [Conventions](#conventions)).

<!-- BEGIN GENERATED: details -->

<a id="details-psi"></a>

### Induction equation (`var_psi`)

Source: `induction_equation_1` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\check{T}_e &= \begin{cases} T_e & \text{if } \texttt{with TiTe} \\ T & \text{otherwise} \end{cases} \\
\hat{T}_e &= \begin{cases} T_e & \text{if } \texttt{with TiTe} \\ \frac{T}{2} & \text{otherwise} \end{cases} \\
p_e &= \rho\,\hat{T}_e + \rho_{imp}\,A_e(\hat{T}_e) \\
w_V &= R\,\mathcal{J}
\end{aligned}
$$

#### Time derivative

$$
A = \frac{1}{R^{2}}\,v\,\psi\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### -B.grad u

$$
\begin{aligned}
&+ \frac{1}{R}\,v\,[\psi,u]^{st} \\
&- \frac{1}{R^{2}}\,v\,F_0\,\partial_{\phi} u
\end{aligned}
$$

##### Eta*j

$$
\begin{aligned}
&+ \frac{1}{R^{2}}\,v\,\eta(\check{T}_e, \rho, \rho_{imp})\,\left(j - j_{src} - j_b\right)
\end{aligned}
$$

##### Hyper-resistivity

$$
\begin{aligned}
&+ \frac{1}{R}\,\eta_{num}(\check{T}_e)\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} j\right)
\end{aligned}
$$

##### Runaway current coupling, -eta*j_RE

$$
\begin{aligned}
&- \frac{1}{R^{2}}\,v\,\eta(\check{T}_e, \rho, \rho_{imp})\,j_{RE}^{ind}
\end{aligned}
$$

##### B.grad Pe diamagnetic term (no impurity contribution)

$$
\begin{aligned}
&- \frac{2}{\rho^{\mathrm{c}}\,B^2(\psi)\,R^{3}}\,v\,\tau_{IC}\,F_0^{2}\,[\psi,p_e]^{st} \\
&+ \frac{2}{\rho^{\mathrm{c}}\,B^2(\psi)\,R^{4}}\,v\,\tau_{IC}\,F_0^{3}\,\partial_{\phi} p_e
\end{aligned}
$$

<a id="details-u"></a>

### Perpendicular momentum equation (`var_u`)

Source: `momentum_equation_2` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\check{T}_e &= \begin{cases} T_e & \text{if } \texttt{with TiTe} \\ T & \text{otherwise} \end{cases} \\
\alpha_e &= \alpha_e(\check{T}_e) \\
S_{ion} &= S_{ion}(\check{T}_e) \\
S_{rec} &= S_{rec}(\check{T}_e) \\
w_V &= R\,\mathcal{J} \\
&\text{if } \texttt{with TiTe}\text{:} \\
\quad p &= \rho\,\left(T_i + T_e\right) + \rho_{imp}\,\left(\alpha_i\,T_i + A_e(T_e)\right) \\
&\text{otherwise:} \\
\quad p &= \rho\,T + \rho_{imp}\,A_{imp}(T) \\
p_i &= p_i^{dia} \\
n_e &= \rho + \alpha_e\,\rho_{imp} \\
S_{neut} &= n_e\,\rho_n\,S_{ion} - n_e\,\left(\rho - \rho_{imp}\right)\,S_{rec}
\end{aligned}
$$

#### Time derivative

$$
A = R^{2}\,w_V\,\left(-\overline{\rho}^{\mathrm{c}}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right) - f_{cons}\,\rho\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \overline{u}\right)\right)
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Perpendicular inertia (poloidal Jacobian of the kinetic energy)

$$
\begin{aligned}
&- \frac{1}{2}\,R\,\left(\left(\partial_{R} u\right)^{2} + \left(\partial_{Z} u\right)^{2}\right)\,[v,\hat{\rho}]
\end{aligned}
$$

##### Vorticity advection

$$
\begin{aligned}
&- R^{3}\,\rho\,\omega\,[v,u]^{st}
\end{aligned}
$$

##### JxB force term. which here is B.gra j

$$
\begin{aligned}
&+ \frac{1}{R}\,v\,[\psi,j]^{st} \\
&- \frac{1}{R^{2}}\,v\,F_0\,\partial_{\phi} j
\end{aligned}
$$

##### Pressure advection

$$
\begin{aligned}
&+ R\,[v,p]^{st}
\end{aligned}
$$

##### Perpendicular and toroidal viscosity

$$
\begin{aligned}
&- \mu_\perp(\check{T}_e)\,R^{2}\,f_\mu^{old}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \omega\right) \\
&- 2\,\mu_\perp(\check{T}_e)\,R\,f_\mu^{new}\,\omega\,\partial_{R} v \\
&- \mu_\perp(\check{T}_e)\,f_\mu^{new}\,\left(\partial_{R} v\,\partial_{R}\partial_{\phi}\partial_{\phi} u + \partial_{Z} v\,\partial_{Z}\partial_{\phi}\partial_{\phi} u\right) \\
&- \frac{1}{R}\,\mu_{num}(\check{T}_e)\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \omega
\end{aligned}
$$

##### Diamagnetic pressure advection and diamagnetic viscosity

$$
\begin{aligned}
&- 2\,v\,\tau_{IC}\,R^{3}\,[p_i,\omega]^{st} \\
&- 2\,\tau_{IC}\,R^{2}\,\partial_{Z} p_i\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right) \\
&- 2\,v\,\tau_{IC}\,R^{3}\,\left(\partial_{R}\partial_{Z} u\,\left(\partial_{R}\partial_{R} p_i - \partial_{Z}\partial_{Z} p_i\right) - \partial_{R}\partial_{Z} p_i\,\left(\partial_{R}\partial_{R} u - \partial_{Z}\partial_{Z} u\right)\right) \\
&+ \mu_\perp'(\check{T}_e)\,W_{dia}\,\left(\nabla_{\mathrm{pol}} \begin{cases} T_i & \text{if } \texttt{with TiTe} \\ \frac{T}{2} & \text{otherwise} \end{cases}\cdot\nabla_{\mathrm{pol}} v\right) \\
&+ \mu_\perp(\check{T}_e)\,W_{dia}\,\nabla^2_{\mathrm{pol}} v
\end{aligned}
$$

##### Conservative form of the momentum equation

$$
\begin{aligned}
&+ f_{cons}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right)\,\left(-R\,[\hat{\rho},u] + F_0\,\left(\rho\,\partial_{\phi} v_\parallel + v_\parallel\,\partial_{\phi} \rho\right) + R\,\rho\,[v_\parallel,\psi] + R\,v_\parallel\,[\rho,\psi]\right)
\end{aligned}
$$

##### Runaway/auxiliary pressure contributions

$$
\begin{aligned}
&- v\,\left(P_\parallel^{RE} + P_\perp^{RE}\right) \\
&+ R\,\left(-\Pi_R^{RE}\,\partial_{Z} v + \Pi_Z^{RE}\,\partial_{R} v\right)
\end{aligned}
$$

##### Momentum carried by the neutral and particle sources

$$
\begin{aligned}
&+ R^{2}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u\right)\,\left(\left(1 - \delta_{n}\right)\,S_{neut} + \left(1 - f_{cons}\right)\,\left(S_p + S_{pellet} + S_{bg,drift} + S_{imp,drift}\right)\right)
\end{aligned}
$$

##### Taylor-Galerkin stabilization

$$
\begin{aligned}
&- \frac{1}{4}\,c_{TG}^{u}\,R^{4}\,\rho\,[\omega,u]\,[v,u]\,\Delta t \\
&- \frac{1}{4}\,c_{TG}^{u}\,R^{2}\,\omega\,f_{cons}\,[\hat{\rho},u]\,[v,u]\,\Delta t
\end{aligned}
$$

**Note:** Neo-classical terms still missing!

<a id="details-zj"></a>

### Current definition (`var_zj`)

Source: `current_constraint_equation_zj` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Constraint

$$
C = \frac{1}{R}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \psi + v\,j\right)\,\mathcal{J}
$$

<a id="details-w"></a>

### Vorticity definition (`var_w`)

Source: `vorticity_constraint_equation_w` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Constraint

$$
C = \left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} u + v\,\omega\right)\,R\,\mathcal{J}
$$

<a id="details-rho"></a>

### Density equation (`var_rho`)

Source: `density_equation_rho` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\check{T}_e &= \begin{cases} T_e & \text{if } \texttt{with TiTe} \\ T & \text{otherwise} \end{cases} \\
\rho_{main} &= \rho - \rho_{imp} \\
n_e &= \rho + \alpha_e(\check{T}_e)\,\rho_{imp} \\
D_\parallel^{\mathrm{tot}} &= D_\parallel + D_{\parallel,sc}\,\tau_{sc} \\
D_{\parallel,imp}^{\mathrm{tot}} &= D_{\parallel,imp} + D_{\parallel,imp,sc}\,\tau_{sc} \\
w_V &= R\,\mathcal{J}
\end{aligned}
$$

#### Time derivative

$$
A = v\,\rho\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Particle sources

$$
\begin{aligned}
&+ v\,\left(S_p + S_{pellet} + S_{bg,drift} + S_{imp,drift} + S_\rho^{aux}\right)
\end{aligned}
$$

##### Sources/sinks due to neutrals

$$
\begin{aligned}
&+ v\,n_e\,\rho_n\,S_{ion}(\check{T}_e) \\
&- v\,n_e\,\rho_{main}\,S_{rec}(\check{T}_e)
\end{aligned}
$$

##### Convection/compression by u variable (ExB)

$$
\begin{aligned}
&+ v\,\mathcal{C}_u(\rho)
\end{aligned}
$$

##### Parallel convection/compression by vpar

$$
\begin{aligned}
&+ v\,\mathcal{C}_\parallel(\rho)
\end{aligned}
$$

##### Parallel and perpendicular diffusion of main ions

$$
\begin{aligned}
&+ \left(D_\parallel^{\mathrm{tot}} - D_\perp\right)\,\mathcal{D}_\parallel(v, \rho_{main}) \\
&- D_\perp\,\mathcal{D}_{tot}(v, \rho_{main})
\end{aligned}
$$

##### Parallel and perpendicular diffusion of impurities (total mass)

$$
\begin{aligned}
&+ \left(D_{\parallel,imp}^{\mathrm{tot}} - D_{\perp,imp}\right)\,\mathcal{D}_\parallel(v, \rho_{imp}) \\
&- D_{\perp,imp}\,\mathcal{D}_{tot}(v, \rho_{imp})
\end{aligned}
$$

##### Numerical stabilization

$$
\begin{aligned}
&+ c_{TG}^{\rho}\,\mathcal{T}_n(v, \rho) \\
&- D_{\perp,num}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \rho
\end{aligned}
$$

##### Diamagnetic drift

$$
\begin{aligned}
&+ 4\,v\,\tau_{IC}\,\partial_{Z} p_i^{dia}
\end{aligned}
$$

##### Pinch term, not sure about this one here

$$
\begin{aligned}
&- \frac{1}{\sqrt{\left(\partial_{R} \psi\right)^{2} + \left(\partial_{Z} \psi\right)^{2}}}\,V_{pinch}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}} \psi\right)\,\rho
\end{aligned}
$$

<a id="details-vpar"></a>

### Parallel velocity equation (`var_vpar`)

The equation is the projection of the momentum equation on $\mathbf{B}$. With the parallel flow $v_\parallel\mathbf{B}$ and $\mathbf{B}\cdot\partial_t\mathbf{B} = \frac{1}{2}\partial_t B^2 = \frac{1}{2}\partial_t B^2_{pol}$, $\mathbf{B}\cdot\rho\,\partial_t(v_\parallel\mathbf{B}) = \rho\,\big(B^2\,\partial_t v_\parallel + \frac{1}{2}v_\parallel\,\partial_t B^2_{pol}\big)$: the time derivative acts separately on $v_\parallel$ and on $B^2_{pol}$, and is not the derivative of $\rho\,v_\parallel B^2$. In the conservative form ($f_{cons} = 1$) it adds $v_\parallel B^2\,\partial_t\rho$. The element routine uses the corrected density $\rho^{\mathrm{c}}$ in the first two terms.

Source: `parallel_velocity_equation_vpar` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\check{T}_e &= \begin{cases} T_e & \text{if } \texttt{with TiTe} \\ T & \text{otherwise} \end{cases} \\
B^2 &= B^2(\psi) \\
w_V &= R\,\mathcal{J} \\
n_e &= \rho + \alpha_e(\check{T}_e)\,\rho_{imp} \\
p &= \begin{cases} \rho\,\left(T_i + T_e\right) & \text{if } \texttt{with TiTe} \\ \rho\,T & \text{otherwise} \end{cases} \\
&\text{if } \texttt{with TiTe}\text{:} \\
\quad p &\mathrel{+}= \rho_{imp}\,\left(\alpha_i\,T_i + A_e(T_e)\right) \\
&\text{otherwise:} \\
\quad p &\mathrel{+}= \rho_{imp}\,A_{imp}(T) \\
\mu_\parallel^{\mathrm{eff}} &= \mu_\parallel + \mu_{\parallel,sc}\,\tau_{sc} \\
S_n &= S_p + S_{pellet} + S_{bg,drift} + S_{imp,drift} \\
\nabla_{\mathrm{pol}}(v_\parallel - V_{rot}) &= \nabla_{\mathrm{pol}} v_\parallel - \frac{\mathrm{d}V_{rot}}{\mathrm{d}\psi}\,\nabla_{\mathrm{pol}} \psi
\end{aligned}
$$

#### Time derivative

$$
A = v\,\left(\overline{\rho}^{\mathrm{c}}\,v_\parallel\,B^2(\overline{\psi}) + \frac{1}{2}\,\overline{\rho}^{\mathrm{c}}\,\overline{v_\parallel}\,B^2_{pol}(\psi) + f_{cons}\,\rho\,\overline{v_\parallel}\,B^2(\overline{\psi})\right)\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Parallel pressure gradient

$$
\begin{aligned}
&- v\,\left(\mathbf{B}\cdot\nabla p\right)
\end{aligned}
$$

##### Parallel advection of the kinetic energy, 0.5*v_par**2*B**2

$$
\begin{aligned}
&+ \frac{1}{2}\,v_\parallel^{2}\,B^2\,\left(\rho\,\left(\mathbf{B}\cdot\nabla v\right) + v\,\left(\mathbf{B}\cdot\nabla \rho\right)\right)
\end{aligned}
$$

##### Term to obtain conservative form of momentum equation

$$
\begin{aligned}
&+ f_{cons}\,v\,v_\parallel\,B^2\,\left(\frac{[\hat{\rho},u]}{R} - \mathbf{B}\cdot\nabla \left(\rho\,v_\parallel\right)\right)
\end{aligned}
$$

##### Numerical and physical parallel viscosities

$$
\begin{aligned}
&- \mu_{\parallel,num}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} v_\parallel \\
&- \frac{1}{R^{2}\,B^2}\,\mu_{\parallel\parallel}\,F_0^{2}\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)\,\left(\mathbf{B}\cdot\nabla v\right) \\
&- \mu_\parallel^{\mathrm{eff}}\,\left(\nabla_{\mathrm{pol}} v\cdot\nabla_{\mathrm{pol}}(v_\parallel - V_{rot})\right)
\end{aligned}
$$

##### External momentum sources

$$
\begin{aligned}
&+ v\,M_\parallel^{aux}
\end{aligned}
$$

##### Momentum carried by the particle sources; not active in conservative form

$$
\begin{aligned}
&- v\,S_n\,v_\parallel\,B^2\,\left(1 - f_{cons}\right) \\
&- v\,S_\rho^{aux}\,v_\parallel\,B^2\,\left(1 - f_{cons}\right)
\end{aligned}
$$

##### This one should probably be multiplied by (1 - fact_conservative_u

$$
\begin{aligned}
&+ \left(1 - \delta_{n}\right)\,v\,n_e\,v_\parallel\,B^2\,\left(\left(\rho - \rho_{imp}\right)\,S_{rec}(\check{T}_e) - \rho_n\,S_{ion}(\check{T}_e)\right)
\end{aligned}
$$

##### Taylor-Galerkin stabilization

$$
\begin{aligned}
&- \frac{1}{4}\,c_{TG}^{v_\parallel}\,\rho\,v_\parallel^{2}\,B^2\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)\,\left(\mathbf{B}\cdot\nabla v\right)\,\Delta t \\
&- \frac{1}{4}\,c_{TG}^{v_\parallel}\,v\,v_\parallel^{2}\,B^2\,\left(1 - f_{cons}\right)\,\left(\mathbf{B}\cdot\nabla v_\parallel\right)\,\left(\mathbf{B}\cdot\nabla \rho\right)\,\Delta t \\
&- \frac{1}{4}\,c_{TG}^{v_\parallel}\,v_\parallel^{3}\,B^2\,f_{cons}\,\left(\mathbf{B}\cdot\nabla \rho\right)\,\left(\mathbf{B}\cdot\nabla v\right)\,\Delta t
\end{aligned}
$$

##### Not sure this pinch term should be here

$$
\begin{aligned}
&+ \frac{1}{\sqrt{\left(\partial_{R} \psi\right)^{2} + \left(\partial_{Z} \psi\right)^{2}}}\,V_{pinch}\,\left(\nabla_{\mathrm{pol}} \psi\cdot\nabla_{\mathrm{pol}} v_\parallel\right)\,\rho\,v
\end{aligned}
$$

<a id="details-rhoimp"></a>

### Impurity density equation (`var_rhoimp`)

`with_TiTe` has no effect here; it is accepted so that every density equation has the same interface.

Source: `impurity_density_equation_rhoimp` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
D_{\parallel,imp}^{\mathrm{tot}} &= D_{\parallel,imp} + D_{\parallel,imp,sc}\,\tau_{sc} \\
w_V &= R\,\mathcal{J}
\end{aligned}
$$

#### Time derivative

$$
A = v\,\rho_{imp}\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Particle sources

$$
\begin{aligned}
&+ v\,S_{imp,drift}
\end{aligned}
$$

##### Convection/compression by u variable (ExB)

$$
\begin{aligned}
&+ v\,\mathcal{C}_u(\rho_{imp})
\end{aligned}
$$

##### Parallel convection/compression by vpar

$$
\begin{aligned}
&+ v\,\mathcal{C}_\parallel(\rho_{imp})
\end{aligned}
$$

##### Parallel and perpendicular diffusion

$$
\begin{aligned}
&+ \left(D_{\parallel,imp}^{\mathrm{tot}} - D_{\perp,imp}\right)\,\mathcal{D}_\parallel(v, \rho_{imp}) \\
&- D_{\perp,imp}\,\mathcal{D}_{tot}(v, \rho_{imp})
\end{aligned}
$$

##### Numerical stabilization

$$
\begin{aligned}
&+ c_{TG}^{\rho_{imp}}\,\mathcal{T}_n(v, \rho_{imp}) \\
&- D_{\perp,num}^{n}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \rho_{imp}
\end{aligned}
$$

<a id="details-rhon"></a>

### Neutral density equation (`var_rhon`)

Fluid neutrals are diffused with an anisotropic diffusivity, convected with the plasma flow when $\delta_n = 1$, ionized and recombined with the same rates as in the density equation, and fed by a prescribed source.

Source: `neutral_density_equation_rhon` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\check{T}_e &= \begin{cases} T_e & \text{if } \texttt{with TiTe} \\ T & \text{otherwise} \end{cases} \\
n_e &= \rho + \alpha_e(\check{T}_e)\,\rho_{imp} \\
\rho_{main} &= \rho - \rho_{imp} \\
w_V &= R\,\mathcal{J}
\end{aligned}
$$

#### Time derivative

$$
A = v\,\rho_n\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Anisotropic diffusion

$$
\begin{aligned}
&- D_{n,R}\,\partial_{R} v\,\partial_{R} \rho_n \\
&- D_{n,Z}\,\partial_{Z} v\,\partial_{Z} \rho_n \\
&- \frac{1}{R^{2}}\,D_{n,\phi}\,\partial_{\phi} v\,\partial_{\phi} \rho_n
\end{aligned}
$$

##### Convection/compression with the plasma flow

$$
\begin{aligned}
&+ \delta_{n}\,v\,\left(\mathcal{C}_u(\rho_n) + \mathcal{C}_\parallel(\rho_n)\right)
\end{aligned}
$$

##### Ionization and recombination with the main/impurity ions

$$
\begin{aligned}
&- v\,n_e\,\rho_n\,S_{ion}(\check{T}_e) \\
&+ v\,n_e\,\rho_{main}\,S_{rec}(\check{T}_e)
\end{aligned}
$$

##### Neutral source

$$
\begin{aligned}
&+ v\,S_n^{drift}
\end{aligned}
$$

##### Numerical stabilization

$$
\begin{aligned}
&- D_{\perp,num}^{n}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} \rho_n
\end{aligned}
$$

<a id="details-Ti"></a>

### Ion energy equation (`var_Ti`)

Source: `ion_energy_equation_Ti` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\mathcal{J}\,[\cdot,\cdot] &= \begin{cases} \mathcal{J}\,[\cdot,\cdot]^{st} & \text{if } \texttt{st form} \\ (a, b) \mapsto [a,b]\,\mathcal{J} & \text{otherwise} \end{cases} \\
w_V &= R\,\mathcal{J} \\
n_i &= \rho + \alpha_i\,\rho_{imp} \\
p_i &= n_i\,T_i
\end{aligned}
$$

#### Time derivative

$$
A = v\,\left(\rho^{\mathrm{c}} + \alpha_i\,\rho_{imp}^{\mathrm{c}}\right)\,T_i\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Convection/compression by u (ExB)

$$
\begin{aligned}
&+ v\,\mathcal{C}_u^{p}(p_i)
\end{aligned}
$$

##### Parallel convection/compression by vpar

$$
\begin{aligned}
&+ v\,\mathcal{C}_\parallel^{p}(p_i)
\end{aligned}
$$

##### Parallel heat conduction

$$
\begin{aligned}
&+ \left(\kappa_{\parallel,i}(T_i) - \kappa_{\perp,i}(\rho)\right)\,\mathcal{D}_\parallel(v, T_i)
\end{aligned}
$$

##### Perpendicular (total) heat conduction

$$
\begin{aligned}
&- \kappa_{\perp,i}(\rho)\,\mathcal{D}_{tot}(v, T_i)
\end{aligned}
$$

##### Heat sources/sinks

$$
\begin{aligned}
&+ v\,\left(H_i + H_i^{aux}\right)
\end{aligned}
$$

##### Energy exchange with electrons

$$
\begin{aligned}
&+ v\,Q_{ie}(T_i, T_e, \rho, \rho_{imp})
\end{aligned}
$$

##### Neutral recombination sink

$$
\begin{aligned}
&- v\,T_i\,\left(\rho^{\mathrm{c}}\right)^{2}\,S_{rec}(T_e)
\end{aligned}
$$

##### Kinetic energy released by particle sources (friction heating)

$$
\begin{aligned}
&+ \frac{1}{2}\,v\,\left(\Gamma - 1\right)\,\left(v_\parallel^{2}\,B^2(\psi) + R^{2}\,\left(\left(\partial_{R} u\right)^{2} + \left(\partial_{Z} u\right)^{2}\right)\right)\,S_{kin}(T_e)
\end{aligned}
$$

##### Parallel viscous heating

$$
\begin{aligned}
&+ Q_{\mu\parallel}(v)
\end{aligned}
$$

##### Perpendicular (u) viscous heating

$$
\begin{aligned}
&+ Q_{\mu\perp}(v, \mu_{heat}(T_e))
\end{aligned}
$$

##### Energy and momentum from the kinetic neutral/impurity model

$$
\begin{aligned}
&+ Q_{kin}(v)
\end{aligned}
$$

##### Numerical stabilization

$$
\begin{aligned}
&+ \mathcal{T}_p(v, p_i; c_{TG}^{T_i}) \\
&- \kappa_{\perp,num}^{i}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} T_i
\end{aligned}
$$

##### Implicit heat source to avoid negative temperatures

$$
\begin{aligned}
&+ \mathcal{H}(v, \epsilon_i^{\mathrm{floor}}(T_i), T_i^{\mathrm{floor}}(T_i), T_{\mathrm{min}})
\end{aligned}
$$

<a id="details-Te"></a>

### Electron energy equation (`var_Te`)

Source: `electron_energy_equation_Te` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\mathcal{J}\,[\cdot,\cdot] &= \begin{cases} \mathcal{J}\,[\cdot,\cdot]^{st} & \text{if } \texttt{st form} \\ (a, b) \mapsto [a,b]\,\mathcal{J} & \text{otherwise} \end{cases} \\
w_V &= R\,\mathcal{J} \\
p_e &= \rho\,T_e + \rho_{imp}\,A_e(T_e) \\
n_e^{\mathrm{c}} &= \rho^{\mathrm{c}} + \alpha_e(T_e)\,\rho_{imp}^{\mathrm{c}} \\
W_{ion} &= E_{ion}(T_e)\,\rho_{imp} + E_{ion}^{bg}\,\left(\rho - \rho_{imp}\right)
\end{aligned}
$$

#### Time derivative

$$
A = v\,\left(\rho^{\mathrm{c}}\,T_e + \rho_{imp}^{\mathrm{c}}\,A_e(T_e) + \left(\Gamma - 1\right)\,W_{ion}\right)\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Convection/compression by u (ExB)

$$
\begin{aligned}
&+ v\,\mathcal{C}_u^{p}(p_e)
\end{aligned}
$$

##### Parallel convection/compression by vpar

$$
\begin{aligned}
&+ v\,\mathcal{C}_\parallel^{p}(p_e)
\end{aligned}
$$

##### Parallel heat conduction

$$
\begin{aligned}
&+ \left(\kappa_{\parallel,e}(T_e) - \kappa_{\perp,e}(\rho)\right)\,\mathcal{D}_\parallel(v, T_e)
\end{aligned}
$$

##### Perpendicular (total) heat conduction

$$
\begin{aligned}
&- \kappa_{\perp,e}(\rho)\,\mathcal{D}_{tot}(v, T_e)
\end{aligned}
$$

##### Heat sources/sinks

$$
\begin{aligned}
&+ v\,\left(H_e + H_e^{aux} + P_{teleport}\right)
\end{aligned}
$$

##### Ionization sink due to neutrals

$$
\begin{aligned}
&- v\,\xi_{ion}\,\left(\rho + \alpha_e(T_e)\,\rho_{imp}\right)\,\rho_n\,S_{ion}(T_e)
\end{aligned}
$$

##### Energy exchange with ions

$$
\begin{aligned}
&+ v\,Q_{ei}(T_i, T_e, \rho, \rho_{imp})
\end{aligned}
$$

##### Ohmic heating

$$
\begin{aligned}
&+ v\,\left(\Gamma - 1\right)\,\eta_{ohm}(T_e, \rho, \rho_{imp})\,\left(\frac{1}{R}\left(j - j_{RE}\right)\right)^{2}
\end{aligned}
$$

##### Neutral (line) radiation

$$
\begin{aligned}
&- v\,n_e^{\mathrm{c}}\,\rho_n^{\mathrm{c}}\,L_{rays}(T_e)
\end{aligned}
$$

##### Background radiation

$$
\begin{aligned}
&- v\,n_e^{\mathrm{c}}\,\left(\rho^{\mathrm{c}} - \rho_{imp}^{\mathrm{c}}\right)\,L_{cont}(T_e) \\
&- v\,n_e^{\mathrm{c}}\,f_{rad}^{bg}(T_e)
\end{aligned}
$$

##### Impurity radiation

$$
\begin{aligned}
&- v\,n_e^{\mathrm{c}}\,\rho_{imp}^{\mathrm{c}}\,L_{rad}(T_e)
\end{aligned}
$$

##### Ionization potential transport

$$
\begin{aligned}
&+ \mathcal{Q}_{ion}(v, W_{ion}, T_e)
\end{aligned}
$$

##### Numerical stabilization

$$
\begin{aligned}
&+ \mathcal{T}_p(v, p_e; c_{TG}^{T_e}) \\
&- \kappa_{\perp,num}^{e}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} T_e
\end{aligned}
$$

##### Implicit heat source to avoid negative temperatures

$$
\begin{aligned}
&+ \mathcal{H}(v, \epsilon_e^{\mathrm{floor}}(T_e), T_e^{\mathrm{floor}}(T_e), T_{\mathrm{min}})
\end{aligned}
$$

<a id="details-T"></a>

### Single-temperature energy equation (`var_T`)

Source: `total_energy_equation_T` in `util/equation_codegen/src/jorek_equations/model600.py`.

#### Definitions

$$
\begin{aligned}
\mathcal{J}\,[\cdot,\cdot] &= \begin{cases} \mathcal{J}\,[\cdot,\cdot]^{st} & \text{if } \texttt{st form} \\ (a, b) \mapsto [a,b]\,\mathcal{J} & \text{otherwise} \end{cases} \\
w_V &= R\,\mathcal{J} \\
p &= \rho\,T + \rho_{imp}\,A_{imp}(T) \\
n_e^{\mathrm{c}} &= \rho^{\mathrm{c}} + \alpha_e(T)\,\rho_{imp}^{\mathrm{c}} \\
W_{ion} &= E_{ion}(T)\,\rho_{imp} + E_{ion}^{bg}\,\left(\rho - \rho_{imp}\right)
\end{aligned}
$$

#### Time derivative

$$
A = v\,\left(\rho^{\mathrm{c}}\,T + \rho_{imp}^{\mathrm{c}}\,A_{imp}(T) + \left(\Gamma - 1\right)\,W_{ion}\right)\,w_V
$$

#### Right-hand side

$$
B = \Big(\textstyle\sum_k b_k\Big)\,w_V
$$

##### Convection/compression by u (ExB)

$$
\begin{aligned}
&+ v\,\mathcal{C}_u^{p}(p)
\end{aligned}
$$

##### Parallel convection/compression by vpar

$$
\begin{aligned}
&+ v\,\mathcal{C}_\parallel^{p}(p)
\end{aligned}
$$

##### Parallel heat conduction

$$
\begin{aligned}
&+ \left(\kappa_\parallel(T) - \kappa_\perp(\rho)\right)\,\mathcal{D}_\parallel(v, T)
\end{aligned}
$$

##### Perpendicular (total) heat conduction

$$
\begin{aligned}
&- \kappa_\perp(\rho)\,\mathcal{D}_{tot}(v, T)
\end{aligned}
$$

##### Heat sources/sinks

$$
\begin{aligned}
&+ v\,\left(H + H^{aux} + P_{teleport}\right)
\end{aligned}
$$

##### Ionization sink due to neutrals

$$
\begin{aligned}
&- v\,\xi_{ion}\,\left(\rho + \alpha_e(T)\,\rho_{imp}\right)\,\rho_n\,S_{ion}(T)
\end{aligned}
$$

##### Recombination sink due to neutrals

$$
\begin{aligned}
&- \frac{1}{2}\,\left(\Gamma - 1\right)\,v\,T\,\left(\rho^{\mathrm{c}}\right)^{2}\,S_{rec}(T)
\end{aligned}
$$

##### Ohmic heating

$$
\begin{aligned}
&+ v\,\left(\Gamma - 1\right)\,\eta_{ohm}(T, \rho, \rho_{imp})\,\left(\frac{1}{R}\left(j - j_{RE}\right)\right)^{2}
\end{aligned}
$$

##### Neutral (line) radiation

$$
\begin{aligned}
&- v\,n_e^{\mathrm{c}}\,\rho_n^{\mathrm{c}}\,L_{rays}(T)
\end{aligned}
$$

##### Background radiation

$$
\begin{aligned}
&- v\,n_e^{\mathrm{c}}\,\left(\rho^{\mathrm{c}} - \rho_{imp}^{\mathrm{c}}\right)\,L_{cont}(T) \\
&- v\,n_e^{\mathrm{c}}\,f_{rad}^{bg}(T)
\end{aligned}
$$

##### Impurity radiation

$$
\begin{aligned}
&- v\,n_e^{\mathrm{c}}\,\rho_{imp}^{\mathrm{c}}\,L_{rad}(T)
\end{aligned}
$$

##### Ionization potential transport

$$
\begin{aligned}
&+ \mathcal{Q}_{ion}(v, W_{ion}, T)
\end{aligned}
$$

##### Kinetic energy released by particle sources (friction heating)

$$
\begin{aligned}
&+ \frac{1}{2}\,v\,\left(\Gamma - 1\right)\,\left(v_\parallel^{2}\,B^2(\psi) + R^{2}\,\left(\left(\partial_{R} u\right)^{2} + \left(\partial_{Z} u\right)^{2}\right)\right)\,S_{kin}(T)
\end{aligned}
$$

##### Parallel viscous heating

$$
\begin{aligned}
&+ Q_{\mu\parallel}(v)
\end{aligned}
$$

##### Perpendicular (u) viscous heating

$$
\begin{aligned}
&+ Q_{\mu\perp}(v, \mu_{heat}(T))
\end{aligned}
$$

##### Energy and momentum from the kinetic neutral/impurity model

$$
\begin{aligned}
&+ Q_{kin}(v)
\end{aligned}
$$

##### Numerical stabilization

$$
\begin{aligned}
&+ \mathcal{T}_p(v, p; c_{TG}^{T}) \\
&- \kappa_{\perp,num}\,\nabla^2_{\mathrm{pol}} v\,\nabla^2_{\mathrm{pol}} T
\end{aligned}
$$

##### Implicit heat source to avoid negative temperatures

$$
\begin{aligned}
&+ \mathcal{H}(v, \epsilon^{\mathrm{floor}}(T), T^{\mathrm{floor}}(T), T_{\mathrm{min}})
\end{aligned}
$$


<!-- END GENERATED: details -->

## Regenerating this page

The equations are edited in `util/equation_codegen/src/jorek_equations/model600.py`
only; do not edit the generated regions of this page (between the
`BEGIN GENERATED` and `END GENERATED` comments), which are overwritten. The
prose outside them is hand-written. How each Python name is written in LaTeX,
and the titles and descriptions of the operators, are in
`util/equation_codegen/src/jorek_equations/model600_notation.py`.

The tools need Python 3 with `sympy`; on the ITER cluster, run
`module load sympy/1.14.0-gfbf-2025b` first. From `util/equation_codegen`:

```bash
python3 examples/model600_docs.py render   # rewrite the generated regions
python3 examples/model600_docs.py check    # change nothing; exit 1 if the page is stale
```

The unit tests and `run_test.sh` run `check`, so a change to `model600.py`
that is not rendered into this page fails the build.
