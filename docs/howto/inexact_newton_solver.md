---
title: "Run with Non-Linear Time-Stepping (Newton)"
nav_order: 7
parent: "Numerics and Stabilization"
grand_parent: "Howto"
layout: default
render_with_liquid: false
---

# Inexact Newton Solver

## Overview

The standard method of setting up and solving the physical equations in JOREK is to linearize the equations with respect to the physical variables once per time step and proceed as described in [JOREK Solver](../numerics/solver.html) to solve the resulting linear system $Ax=b$. This method may not be valid for large time steps or in the nonlinear phase. It is therefore desirable to have a method that, at least in principle, can solve the physical equations nonlinearly. This can be achieved using an inexact Newton method, which has previously been implemented to some extent in JOREK [reference 1](#ref-1). The solver's sensitivity to whether the current phase of the run is linear or nonlinear also makes it possible to implement adaptive time stepping.

## The Inexact Newton Method

### Newton's method

Write the nonlinear system of time-linearized physical equations as

$$
F(u^{n+1})=F(u^{n+1};u^n,u^{n-1})
=G(u^{n+1})-b(u^n,u^{n-1})=0.
\tag{1}
$$

By abuse of notation, $F$ implicitly depends on $u^n,u^{n-1}$, while $G$ and $b$ depend on the choice of time-integration method through the parameters $\theta$ and $\xi$; see [Time Integration](../numerics/time-integration.html) and page 23 of [reference 1](#ref-1). The standard method linearizes Eq. (1) at $u^n$ and solves

$$
F^{\prime}(u^{n})\delta u^{n+1}
=-F(u^n)
=-G(u^{n})+b(u^n,u^{n-1})=0.
\tag{2}
$$

Here, $\delta u^{n+1}:=u^{n+1}-u^n$ corresponds to the sought solution vector $x$, and

$$
F^{\prime}(u^n)
=\frac{\partial F(u^{n+1};u^n,u^{n-1})}{\partial u^{n+1}}
\bigg\vert_{u^{n+1}=u^n}
=\frac{\partial G(u^{n+1})}{\partial u^{n+1}}
\bigg\vert_{u^{n+1}=u^n}.
\tag{3}
$$

Thus, the standard method linearizes once and solves for $\delta u^{n+1}$ to compute the physical variables at step $n+1$.

Newton's method solves Eq. (1) iteratively. For a given time step $n\to n+1$, compute a series $\{u_k\}_{k\in\mathbb{N}}$ such that $u_k\to u^{n+1}$ as $k\to\infty$. Linearize Eq. (1) at $u_{k=0}=u^n+\delta u_0$, using an initial estimate $\delta u_0$:

$$
F(u_{1})\simeq F(u_0)+(u_{1}-u_0)F^{\prime}(u_0)
\stackrel{!}{=}0
\Rightarrow
F^{\prime}(u_0)\delta u_{1}
=-F(u_0)
=-G(u^{n}+\delta u_0)+b(u^n,u^{n-1}).
\tag{4}
$$

Here, $u_{k=1}=u_{k=0}+\delta u_{k=1}$ is the next estimate for the physical variables at $n+1$. At each iteration $k$, the linear system is

$$
F^{\prime}(u_k)\delta u_{k+1}=-F(u_k),
\tag{5}
$$

where $F(u_k)=F(u_k;u^n,u^{n+1})$.

The Jacobian $F^{\prime}(u^n)=G^{\prime}(u^n)$ and the right-hand side of Eq. (2) are constructed in `call construct_matrix(...)`. In Eq. (4), only $G(u_k)$ needs to be updated. Because $G$ is not accessible independently in the code, write the right-hand side of Eq. (5) as

$$
-F(u_k)
=-G(u_k)+b(u^n,u^{n-1})-G(u^n)+G(u^n)
=-F(u^n)-F^{\prime}(u_k)(u_k-u^n),
\tag{6}
$$

and hence

$$
F^{\prime}(u_k)\delta u_{k+1}
=-F(u^n)-F^{\prime}(u_k)\delta u_k^n,
\qquad
\delta u_k^n=u_k-u^n
=\delta u_k+\delta u_{k-1}+\dots+\delta u_0.
\tag{7}
$$

For each step $k$, the Jacobian $F^{\prime}(u_k)$ is computed through `construct_matrix`. The $k$-th right-hand side is obtained by subtracting $F^{\prime}(u_k)\delta u_k^n$ from the $n$-th right-hand side, after which a solver is called. Once convergence is achieved, $u^{n+1}=u_k=u^n+\delta u_k^n$.

### The “inexact” part

Only an inexact Newton method can be used because the linear system $Ax=b$ in Eq. (7) can only be solved to a finite tolerance. For convergence, the sequence of residuals of the solutions of Eq. (7) must satisfy

$$
\frac{\left|F^{\prime}(u_k)\delta u_{k+1}+F(u_k)\right|}
{\left|F(u_k)\right|}
<1
\quad\Leftrightarrow\quad
\epsilon_k<\left|F(u_k)\right|.
\tag{8}
$$

The tolerances $\{\epsilon_k\}$ must be chosen accordingly; see [reference 2](#ref-2). Otherwise, $\delta u_{k+1}=0$ for all $k$ would be an acceptable solution. Analogous to [reference 1](#ref-1), set

$$
\epsilon_k
=\min\left(\epsilon_0,\frac{1}{2}\epsilon_{k-1},\eta_k\right),
\qquad
\eta_k
:=\gamma\left(\frac{|F(u_k)|}{|F(u_{k-1})|}\right)^\alpha,
\tag{9}
$$

where the usual values are $\gamma=0.5$ and $\alpha=2$. Iterations are only sensible while Eq. (8) is satisfied, or while the right-hand side in Eq. (7) decreases, $|F(u_k)|<|F(u_{k-1})|$. In practice, the computation is aborted if $|F(u_k)|>|F(u_{k-1})|$ for two consecutive iterates.

## Usage

To use the implemented inexact Newton method, check out the relevant feature branch:

```bash
git checkout feature/nlin_tstep
```

Add the following line to `Makefile.inc`:

```makefile
USE_NEWTON = 1
```

The module in `solvers/mod_inexact_newton.f90` is loaded in the main program:

```fortran
program JOREK2
...
#ifdef USE_NEWTON
use mod_inexact_newton, only: inexact_newton
#endif
implicit none
...
```

The `inexact_newton` subroutine is called in the solver section of `jorek2_main.f90`:

```fortran
#ifdef USE_NEWTON
      call inexact_newton(...)
#elif USE_BICGSTAB
      call bicgstab_driver(...)
#else
      call gmres_driver(...)
#endif
```

The behavior of the inexact Newton method is governed by input variables:

```fortran
newton_start       = 5.d-1   !< start value, delta_k = deltas*newton_start
newton_gamma       = 5.d-1   !< gamma for computation of eps_k
newton_alpha       = 2.0d0   !< alpha for computation of eps_k
newton_eps_0       = 1.d-4   !< tolerance for initial Newton step
newton_eps_f       = 1.d-6   !< tolerance for final Newton step
newton_max_iter    = 20      !< maximum number of Newton iterations
```

The upper limit of GMRES iterations is set by `gmres_max_iter` as usual. The variables `newton_alpha` and `newton_gamma` correspond to Eq. (9) and can be left unchanged. The starting value for the $k=0$ step in Eq. (4) is $\delta u_0=\delta u^n\times\mathtt{newton\_start}$.

In most tested cases, `newton_eps_0/newton_eps_f = 10..1000` leads to fast convergence. In general, `newton_eps_f > gmres_tol`, while `newton_eps_f/gmres_tol = 10..100` leads to physically sensible results. Here, `gmres_tol = 1.d-5...1.d-8` is the tolerance that would be set for the standard case. The parameter `newton_max_iter` sets the maximum number of Newton iterations; the run stops once it is reached.

## Adaptive Time Stepping

To use adaptive time stepping, check out the corresponding feature branch:

```bash
git checkout feature/nlin_tstep_adaptive
```

Again, add this line to `Makefile.inc`:

```makefile
USE_NEWTON = 1
```

The additional input parameters are:

```fortran
newton_adapt_time = .false. !< switch for adaptive time stepping
newton_alpha_dt   = 1.05    !< fast convergence: increase tstep
newton_beta_dt    = 0.95    !< slow convergence: decrease tstep
newton_gamma_dt   = 0.75    !< non-convergence: retry with a smaller tstep
```

Activate the feature with:

```fortran
newton_adapt_time = .true.
```

The time step is adapted according to convergence:

- **Fast convergence:** If $\sum_k \mathtt{iter\_gmres}\leq\mathtt{gmres\_max\_iter}/2$ , increase the time step for $n+1\to n+2$ with `tstep = tstep*newton_alpha_dt`.
- **Slow convergence:** If $\sum_k \mathtt{iter\_gmres}>\mathtt{gmres\_max\_iter}/2$, decrease the time step for $n+1\to n+2$ with `tstep = tstep*newton_beta_dt`.
- **No convergence:** Recompute the step $n\to n+1$ with `tstep = tstep*newton_gamma_dt`.

Here, $k$ is the Newton-iteration index. There is a hard-coded minimum time step, `t_min = 0.001`. The computation is aborted if convergence cannot be achieved with `tstep = t_min`. The parameters `newton_alpha_dt`, `newton_beta_dt`, and `newton_gamma_dt` should be actively fine-tuned for the case.

## Limitations

Only `jorek_model = 199, 303, 501, 502, 600, 710, 711, 712` have been included and tested. In principle, a new model can be added by including all new parameters in the corresponding `models/modelXYZ/initialise_parameters.f90`.

## References

1. <a id="ref-1"></a>Franck, E., Hoelzl, M., Lessig, A., & Sonnendrucker, E. (2014). “Energy Conservation and numerical stability for the reduced MHD models of the non-linear JOREK code.” *Mathematical Modelling and Numerical Analysis*, 49, 1331–1365. [arXiv:1408.2099v3](https://arxiv.org/abs/1408.2099)
2. <a id="ref-2"></a>Dembo, R. S., Eisenstat, S. C., & Steihaug, T. (1982). “Inexact Newton Methods.” *SIAM Journal on Numerical Analysis*, 19(2), 400–408. [doi:10.1137/0719025](https://doi.org/10.1137/0719025)
