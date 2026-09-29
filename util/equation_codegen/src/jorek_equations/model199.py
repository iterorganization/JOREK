"""Canonical fields and weak equations for JOREK model 199."""

import sympy as sp

from .equations import ConstraintEquation, EvolutionEquation
from .external import external_function
from .operators import R, bracket, dR, dZ, dphi, dot, laplacian, poiss_bracket_st, grad
from .symbols import coefficient, field, freeze, test_function


psi = field("psi", fortran_current="ps0", fortran_trial="psi")
u = field("u", fortran_current="u0", fortran_trial="u")
j = field("j", fortran_current="zj0", fortran_trial="zj")
omega = field("omega", fortran_current="w0", fortran_trial="w")
rho = field("rho", fortran_current="r0", fortran_trial="rho")
T = field("T", fortran_current="T0", fortran_trial="T")

FIELDS = (psi, u, j, omega, rho, T)


# The constitutive interface is deliberately explicit: the value is evaluated
# by the host application and the tangent is named for generated Fortran.
eta = external_function(
    "eta",
    arguments=("T",),
    derivatives={"T": "deta_dT"},
    policy="piecewise_active",
    fortran_name="eta_T",
)

xjac = coefficient("xjac")
eps_cyl = coefficient("eps_cyl")
F0 = coefficient("F0")
eta_num = coefficient("eta_num")
current_source = coefficient("current_source")
GAMMA = coefficient("GAMMA")
gamma = coefficient("gamma")
visco_num = coefficient("visco_num")
D_par = coefficient("D_par")
ZK_par = coefficient("ZK_par")
D_prof = coefficient("D_prof")
ZK_prof = coefficient("ZK_prof")
particle_source = coefficient("particle_source")
heat_source = coefficient("heat_source")
visco = external_function(
    "visco", arguments=("T",), derivatives={"T": "dvisco_dT"},
    policy="piecewise_active", fortran_name="visco_T"
)
eta_ohmic = external_function(
    "eta_ohmic", arguments=("T",), derivatives={"T": "deta_dT_ohm"},
    policy="piecewise_active", fortran_name="eta_T_ohm"
)


def _rho_hat():
    return R**2 * rho


def _rho_hat_time():
    # The momentum time coefficient is frozen in the Fortran assembly.
    return R**2 * freeze(rho)


def _pressure():
    return rho * T


def _vv2():
    """Squared poloidal velocity used by the momentum equation."""

    return R**2 * (dR(u)**2 + dZ(u)**2)


def _parallel_gradient(value, flux):
    """Full parallel-gradient factor in the physical-coordinate form."""

    return bracket(value, flux) / R + F0 * dphi(value) / R**2


def _parallel_norm():
    return (F0**2 + dR(psi)**2 + dZ(psi)**2) / R**2


def induction_equation_1():
    """Return model-199 equation 1 in its already-integrated weak form.

    The expression is written in physical (R, Z, phi) derivatives.  The
    factor ``xjac`` is explicit because the existing element routine performs
    the poloidal mapping and Gaussian weight outside the equation terms.
    """

    v = test_function("v")
    B = (
        v * eta(T) * (j - current_source) / R * xjac
        + v * xjac * bracket(psi, u)
        - v * eps_cyl * F0 / R * dphi(u) * xjac 
        + eta_num * xjac * dot(grad(v), grad(j))
    )
    A = v * psi / R * xjac
    return EvolutionEquation("model199_induction", v, A, B)


def momentum_equation_2():
    """Return model-199 equation 2 in its integrated weak form."""

    v = test_function("v")
    rho_hat = _rho_hat()
    rho_hat_time = _rho_hat_time()
    pressure = _pressure()
    B = (
        -sp.Rational(1, 2) * _vv2()
        * (dR(v) * dZ(rho_hat) - dZ(v) * dR(rho_hat)) * xjac
        - rho_hat * R**2 * omega * poiss_bracket_st(v, u)
        + v * poiss_bracket_st(psi, j)
        - visco(T) * R * dot(grad(v), grad(omega)) * xjac
        - v * eps_cyl * F0 / R * dphi(j) * xjac
        + R**2 * poiss_bracket_st(v, pressure)
        - visco_num * laplacian(v) * laplacian(omega) * xjac
    )
    A = -R * rho_hat_time * dot(grad(v), grad(u)) * xjac
    return EvolutionEquation("model199_momentum", v, A, B)


def current_constraint_equation_3():
    """Return model-199 equation 3, a static current constraint."""

    v = test_function("v")
    C = (dot(grad(v), grad(psi)) + v * j) / R * xjac
    return ConstraintEquation("model199_current_constraint", v, C, kind="static")


def vorticity_constraint_equation_4():
    """Return model-199 equation 4, a static vorticity constraint."""

    v = test_function("v")
    C = (dot(grad(v), grad(u)) + v * omega) * R * xjac
    return ConstraintEquation("model199_vorticity_constraint", v, C, kind="static")


def density_equation_5():
    """Return model-199 equation 5 (density transport)."""

    v = test_function("v")
    rho_hat = R**2 * rho
    B = (
        v * R * particle_source * xjac
        + v * R**2 * poiss_bracket_st(rho, u)
        + 2 * v * R * rho * dZ(u) * xjac
        - (D_par - D_prof) * R / _parallel_norm()
        * _parallel_gradient(v, psi) * _parallel_gradient(rho, psi) * xjac
        - D_prof * R * dot(grad(v), grad(rho)) * xjac
        - D_prof * R * eps_cyl**2 / R**2 * dphi(v) * dphi(rho) * xjac
    )
    A = v * R * rho * xjac
    return EvolutionEquation("model199_density", v, A, B)


def temperature_equation_6():
    """Return model-199 equation 6 (temperature transport)."""

    v = test_function("v")
    B = (
        v * R * heat_source * xjac
        + v * R**2 * poiss_bracket_st(T, u)
        + 2 * (GAMMA - 1) * v * R * T * dZ(u) * xjac
        - (ZK_par - ZK_prof) * R / _parallel_norm()
        * _parallel_gradient(v, psi) * _parallel_gradient(T, psi) * xjac
        - ZK_prof * R * dot(grad(v), grad(T)) * xjac
        - ZK_prof * R / R**2 * dphi(v) * dphi(T) * xjac
        + v * (gamma - 1) * eta_ohmic(T) * (j / R)**2 * R * xjac
    )
    A = v * R * T * xjac
    return EvolutionEquation("model199_temperature", v, A, B)
