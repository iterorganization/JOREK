# JOREK equation-codegen specification

## 1. Purpose and scope

This document is the initial correctness reference for the Python equation
code generator. The first target is an exact, auditable reproduction of the
volume terms assembled by `models/model199/mod_elt_matrix_fft.f90`.

Version 1 accepts equations that are already in weak form. It does not derive
the weak form from a strong equation and does not perform integration by
parts. Strong-form input, automatic integration by parts, and explicit surface
operators are deferred to a later version.

There are two distinct goals:

1. **Legacy reproduction mode:** reproduce model 199, including quantities
   that the existing Jacobian intentionally or historically treats as frozen.
2. **Declared-physics mode:** later allow each external or derived quantity to
   state its true dependencies explicitly and differentiate them when desired.

The generator must never change from one policy to the other implicitly.

## 2. Model-199 fields

The ordered state vector is

$$
\mathbf q=(\psi,u,j,\omega,\rho,T).
$$

The ordering is part of the generated-code interface:

| Index | DSL name | Mathematical meaning | Current-state Fortran name | Trial-function Fortran name |
|---:|---|---|---|---|
| 1 | `psi` | poloidal magnetic flux | `ps0` | `psi` |
| 2 | `u` | velocity stream function | `u0` | `u` |
| 3 | `j` | toroidal plasma current | `zj0` | `zj` |
| 4 | `omega` | toroidal vorticity | `w0` | `w` |
| 5 | `rho` | particle density | `r0` | `rho` |
| 6 | `T` | total temperature | `T0` | `T` |

The DSL uses `j` and `omega` as the canonical physics names. The Fortran
backend maps them to the legacy names `zj` and `w`.

## 3. State, increment, and trial-function naming

For a field `q`, the DSL distinguishes the following roles:

| DSL concept | Mathematical notation | Example for `psi` | Purpose |
|---|---|---|---|
| field declaration | $q$ | `psi = field("psi")` | identifies a state variable |
| current/background value | $q_0\equiv q^n$ | `psi.current` | evaluates the residual |
| symbolic increment | $\delta q$ | `delta(psi)` | directional linearization |
| previous increment | $\delta q^{n-1}$ | `previous_delta(psi)` | multistep time term |
| finite-element trial function | $\varphi_j$ | `trial(psi)` | matrix column after differentiation |

Linearization is defined by the Gateaux derivative

$$
D F(\mathbf q_0)[\delta\mathbf q]
=
\left.\frac{d}{d\epsilon}
F(\mathbf q_0+\epsilon\,\delta\mathbf q)
\right|_{\epsilon=0}.
$$

When a Jacobian block for field $q$ is emitted, `delta(q)` and all its
derivatives are replaced by the corresponding finite-element trial function
and its derivatives. Background values keep the legacy Fortran suffix `0`.

Equation definitions use the bare field declaration (`psi`, `u`, `rho`, ...),
which denotes an abstract state field. The `.current` form is an explicit
low-level state value and is useful in tests or when constructing already
linearized expressions. The linearization stage resolves every non-varied
abstract field to its current/background value and replaces only the selected
field by the requested increment or trial role.

Spatial variation and differentiation commute:

$$
\delta(\partial_a q)=\partial_a(\delta q),
\qquad a\in\{R,Z,\phi\}.
$$

## 4. Coordinates and geometry

### 4.1 Physical coordinates

The physical coordinates are cylindrical

$$
(R,Z,\phi).
$$

The DSL coordinate names are `R`, `Z`, and `phi`, with derivative operators
`dR`, `dZ`, and `dphi`. The legacy model-199 implementation uses `x`, `y`, and
`p` suffixes for these derivatives:

| DSL | Mathematical notation | Legacy suffix |
|---|---|---|
| `dR(f)` | $\partial_R f$ | `f_x` |
| `dZ(f)` | $\partial_Z f$ | `f_y` |
| `dphi(f)` | $\partial_\phi f$ | `f_p` |

The toroidal direction is periodic. Toroidal endpoint surface terms therefore
cancel.

### 4.2 Element coordinates

The two-dimensional element coordinates are $(s,t)$, with mapping

$$
R=R(s,t),\qquad Z=Z(s,t),
$$

and poloidal Jacobian

$$
J=R_s Z_t-R_t Z_s.
$$

The physical derivatives are

$$
\partial_R f=\frac{Z_t f_s-Z_s f_t}{J},\qquad
\partial_Z f=\frac{-R_t f_s+R_s f_t}{J}.
$$

The DSL's canonical equations are written using physical derivatives. The
Fortran backend may rewrite antisymmetric products in $(s,t)$ coordinates,
as the existing routine does, but that is a code-generation transformation and
must not alter the symbolic equation.

The weak integrals use $dR\,dZ\,d\phi$. Cylindrical factors of $R$ remain
explicit in the equation. At a quadrature point, $dR\,dZ$ contributes the
mapping factor `xjac`, and Gaussian integration contributes `wst`.

Geometry is frozen during physics linearization: variations of `R`, `Z`, the
mapping, `xjac`, metric factors, basis functions, and quadrature weights are
zero.

## 5. Test- and trial-function conventions

The default discretization is Galerkin: test and trial functions come from the
same finite-element space but have independent row and column indices.

```python
v = test_function("v")
```

denotes a symbolic test function. It is not a physical field and is never
linearized. Its derivatives `dR(v)`, `dZ(v)`, and `dphi(v)` are valid symbolic
objects.

For a residual component $\mathcal R_i$, the local row uses test basis
function $v_i$. A variation is expanded as

$$
\delta q=\sum_j \varphi_j\,\delta q_j,
$$

and the matrix entry is

$$
M_{ij}^{(q)}=D\mathcal R_i[\varphi_j].
$$

This corresponds to the outer `(i,j)` basis loops for `v` and the inner
`(k,l)` basis loops for the trial quantities in the current Fortran routine.

Toroidal derivatives determine the intermediate FFT channel:

| Derivative placement | RHS channel | Matrix channel |
|---|---|---|
| neither test nor trial | `RHS_p` | `ELM_p` |
| trial only | not applicable | `ELM_n` |
| test only | `RHS_k` | `ELM_k` |
| both test and trial | not applicable | `ELM_kn` |

The names `p`, `n`, `k`, and `kn` are legacy output names. Internally, the DSL
must retain semantic flags such as `test_phi_order` and `trial_phi_order`.
Unsupported derivative orders must produce an error instead of being assigned
heuristically.

## 6. Equation, RHS, and AMAT sign conventions

### 6.1 Evolution equations

Evolution equations 1, 2, 5, and 6 are represented as weak functionals

$$
\frac{\partial A_i(\mathbf q)}{\partial t}=B_i(\mathbf q).
$$

The state update is

$$
\mathbf q^{n+1}=\mathbf q^n+\delta\mathbf q^n.
$$

With JOREK's `theta` and `zeta` convention, the generated local system is

$$
\left[(1+\zeta)D A_i
-\theta\,\Delta t\,D B_i\right][\delta\mathbf q^n]
=
\Delta t\,B_i(\mathbf q^n)
+\zeta\,D A_i[\delta\mathbf q^{n-1}].
$$

Therefore:

```text
RHS  = + tstep * B(current)
       + zeta * D(A)(previous_delta)

AMAT = + (1 + zeta) * D(A)(trial)
       - theta * tstep * D(B)(trial)
```

The existing routine rescales the stored previous increment by
`tstep/tstep_prev` before using it. This rescaling belongs to time-history
preparation, not symbolic differentiation.

For exact legacy reproduction, the momentum time functional treats the
current-state density coefficient `r0_hat` as frozen while differentiating the
time term with respect to `u`. This special policy must be declared explicitly;
it must not arise from generic differentiation.

### 6.2 Algebraic constraint equations

Equations 3 and 4 are algebraic definitions rather than evolution equations.
For a weak constraint

$$
C_i(\mathbf q)=0,
$$

the Newton correction convention is

$$
D C_i(\mathbf q^n)[\delta\mathbf q]=-C_i(\mathbf q^n).
$$

Therefore:

```text
RHS  = -C(current)
AMAT = +D(C)(trial)
```

No `theta`, `zeta`, or `tstep` factor is applied to these constraint blocks.

## 7. Boundary-term policy

Version 1 starts from the already-integrated weak form. Its input expression is
the volume integrand that remains after any integration by parts. It does not
reconstruct a discarded surface term.

For legacy model-199 reproduction:

1. Toroidal surface terms cancel because $\phi$ is periodic.
2. Poloidal surface terms produced by integrations by parts are absent from
   `mod_elt_matrix_fft.f90` and are treated as discarded in the volume
   generator.
3. The model applies essential boundary conditions separately after volume
   assembly. Fixed-boundary test/variation traces are consequently constrained.
4. Free-boundary or open-boundary surface physics is outside the version-1
   generator scope. The absence of a surface contribution must not be taken as
   a derivation that it is zero for those cases.

The discarded poloidal terms belong to these classes:

| Weak volume structure | Discarded surface-term class |
|---|---|
| `grad(v) . grad(psi)` and `grad(v) . grad(u)` | normal flux from the elliptic current/vorticity definitions |
| `grad(v) . grad(j)` and `grad(v) . grad(omega)` | normal numerical resistive/viscous flux |
| `grad(v) . grad(rho)` and `grad(v) . grad(T)` | normal perpendicular particle/heat flux |
| parallel-gradient test factor times parallel field gradient | normal anisotropic particle/heat flux |
| derivatives moved onto `v` in bracket or pressure terms | advective or pressure boundary flux generated by that integration by parts |
| the fourth-order numerical-viscosity weak product | boundary traces generated by the repeated integration by parts |

Because version 1 receives only the final weak expression, the exact signed
surface formula is not inferred. Each future strong-form transformation must
return both a volume term and an explicit boundary term. Discarding that term
will require a named boundary policy and a stated boundary condition.

## 8. Linearization dependency policy

Every nontrivial derived quantity or external function has one of three
policies:

- `active`: apply the chain rule using its declared derivatives;
- `frozen`: evaluate its value at the current state and set its variation to
  zero;
- `piecewise_active`: select a branch at the current state, then use derivatives
  explicitly supplied for that branch.

If an active dependency has no derivative interface, generation must fail with
a diagnostic. It must not silently freeze the dependency.

### 8.1 Quantities differentiated in legacy reproduction mode

| Quantity | Active dependencies | Required derivative behavior |
|---|---|---|
| all six state fields and their spatial derivatives | their own field | standard linear variation |
| $P=\rho T$ | `rho`, `T` | product rule |
| $\hat\rho=R^2\rho$ and its gradients | `rho` | geometry remains frozen |
| $\lvert\mathbf v_E\rvert^2=R^2\lvert\nabla u\rvert^2$ (`vv2`) | `u` | differentiate both gradient factors |
| $B^2=(F_0^2+\lvert\nabla\psi\rvert^2)/R^2$ | `psi` | includes the derivative called `BB2_psi` |
| parallel-gradient expressions | `psi` and the transported field | product and quotient rules, including the derivative of `B^2` |
| `eta_T` | `T` | external interface supplies `deta_dT` |
| `eta_T_ohm` | `T` | external interface supplies `deta_dT_ohm` |
| `visco_T` | `T` | external interface supplies `dvisco_dT` |

Temperature-dependent coefficients use `piecewise_active`: the current branch
or clipping regime is selected first, and the interface supplies the derivative
for that branch. A clipped constant branch consequently has zero derivative.

### 8.2 Quantities frozen in legacy reproduction mode

| Quantity | Notes |
|---|---|
| coordinates, mapping, metric, `xjac`, basis functions, and quadrature weights | mesh/shape differentiation is out of scope |
| `F0`, `gamma`, `theta`, `zeta`, `tstep`, and numerical coefficients | run parameters |
| `D_par`, `ZK_par`, `eta_num`, and `visco_num` | scalar coefficients |
| `psi_norm = get_psi_n(psi, Z)` | evaluated from the current state but not differentiated in model 199 |
| `D_prof = get_dperp(psi_norm)` | frozen even though its evaluated value depends on `psi_norm` |
| `ZK_prof = get_zkperp(psi_norm)` | frozen even though its evaluated value depends on `psi_norm` |
| threshold-based replacements `D_prof_neg` and `ZK_prof_neg` | branch and value frozen during a Jacobian evaluation |
| current, particle, and heat sources | frozen, including any state dependence inside their evaluation routines |
| equilibrium and normalization data such as axis/boundary fluxes | external data |

Model 199 evaluates the background density and temperature with `abs(...)`,
but its matrix expressions use `rho` and `T` without the derivative of the
absolute-value operation. To reproduce the implementation, background
sanitization is treated as an external evaluation step while field increments
use the identity tangent. This behavior must be marked as a legacy policy.

### 8.3 Future explicit external-function interface

An external function declaration must record:

```text
name
arguments and declared state dependencies
value symbol or value callback
available first derivatives
linearization policy
piecewise/clipping policy, when applicable
Fortran names for the value and derivatives
```

For example, the mathematical function `eta(T)` supplies both its value and
`deta_dT`. Profile functions can later be changed from `frozen` to `active`
only by an explicit specification change and the addition of the corresponding
derivative interface.

## 9. Correctness requirements for later implementation

The eventual Python implementation must:

1. preserve field, test, trial, and coordinate-derivative roles in the symbolic
   expression tree;
2. reproduce the sign rules in Section 6 before applying algebraic
   simplification;
3. make every external dependency policy inspectable;
4. keep named physics terms available in generated comments and audit output;
5. classify toroidal derivatives structurally;
6. reject missing derivatives, unsupported toroidal orders, and undeclared
   boundary-term deletion;
7. compare generated model-199 terms against the existing RHS and AMAT
   expressions using symbolic and numerical tests.

## 10. Scope decisions still requiring physics review

The following legacy behaviors are recorded for reproduction but should be
reviewed before declaring them the preferred mathematical model:

- freezing `D_prof`, `ZK_prof`, `psi_norm`, and state-dependent sources;
- freezing the density coefficient in the momentum time derivative;
- using an identity tangent after `abs(...)` background sanitization;
- omitting poloidal surface contributions outside fixed-boundary use.

Changes to these items define a different Jacobian or boundary-value problem
and must be made as explicit specification revisions.
