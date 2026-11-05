"""Model-600 weak equations, extending the model-199 DSL definitions."""

import functools

import sympy as sp

from .equations import ConstraintEquation, EvolutionEquation
from .external import external_function
from .model199 import FIELDS as MODEL199_FIELDS
from .model199 import GAMMA, gamma, j, omega, psi, rho, T, u, xjac
from .operators import R, dR, dZ, dphi, ds, dt, poiss_bracket_st, poiss_bracket, laplacian, dot, grad
from .symbols import coefficient, field, freeze, test_function


# Model-600 extension fields.  The model199 declarations are reused so their
# current/trial and Fortran naming conventions remain identical.
# -------------------------------------------------------------------------
# Fields
# -------------------------------------------------------------------------

vpar   = field("vpar", fortran_current="vpar0", fortran_trial="vpar")
Ti     = field("Ti", fortran_current="Ti0", fortran_trial="Ti")
Te     = field("Te", fortran_current="Te0", fortran_trial="Te")
rhon   = field("rhon", fortran_current="rn0", fortran_trial="rhon")
rhoimp = field("rhoimp", fortran_current="rimp0", fortran_trial="rhoimp")

FIELDS = MODEL199_FIELDS + (vpar, Ti, Te, rhon, rhoimp)

# -------------------------------------------------------------------------
# Scalar work values
# -------------------------------------------------------------------------
# Values the element routine computes once per quadrature point and
# freezes in the Newton tangent.
Jb = coefficient("Jb")
aux_jre_ind = coefficient("aux_jre_ind")
tauIC = coefficient("tauIC")
tstep = coefficient("tstep")
r0_corr = coefficient("r0_corr")
F0 = coefficient("F0")
factor_psi = tuple(coefficient("factor_psi_{}".format(index)) for index in range(1, 7))
visco_fact_old = coefficient("visco_fact_old")
visco_fact_new = coefficient("visco_fact_new")
dvisco_dT = coefficient("dvisco_dT")
fact_conservative_u = coefficient("fact_conservative_u")
delta_n_convection = coefficient("delta_n_convection")
tgnum_u = coefficient("tgnum_u")
amu_neo_prof = coefficient("amu_neo_prof")
aki_neo_prof = coefficient("aki_neo_prof")
Btheta2 = coefficient("Btheta2")
epsil = coefficient("epsil")
alpha_e = coefficient("alpha_e")
alpha_i = coefficient("alpha_i")
alpha_imp = coefficient("alpha_imp")
Sion_T = coefficient("Sion_T")
Srec_T = coefficient("Srec_T")
particle_source = coefficient("particle_source")
source_pellet = coefficient("source_pellet")
source_bg_drift = coefficient("source_bg_drift")
source_imp_drift = coefficient("source_imp_drift")
aux_P_par_re = coefficient("aux_P_par_re")
aux_P_perp_re = coefficient("aux_P_perp_re")
aux_divPIR_perp = coefficient("aux_divPIR_perp")
aux_divPIZ_perp = coefficient("aux_divPIZ_perp")
D_par_local = coefficient("D_par_local")
D_par_local_imp = coefficient("D_par_local_imp")
D_par_sc_num = coefficient("D_par_sc_num")
D_par_imp_sc_num = coefficient("D_par_imp_sc_num")
D_perp_num_psin = coefficient("D_perp_num_psin")
D_prof = coefficient("D_prof")
D_prof_imp = coefficient("D_prof_imp")
tau_sc = coefficient("tau_sc")
tgnum_rho = coefficient("tgnum_rho")
V_prof_pinch = coefficient("V_prof_pinch")
aux_rho0 = coefficient("aux_rho0")
visco_par = coefficient("visco_par")
visco_par_sc_num = coefficient("visco_par_sc_num")
visco_par_num = coefficient("visco_par_num")
visco_par_par = coefficient("visco_par_par")
tgnum_vpar = coefficient("tgnum_vpar")
aux_mom_par0 = coefficient("aux_mom_par0")
dV_dpsi_source = coefficient("dV_dpsi_source")
Dn_perp_num = coefficient("Dn_perp_num")
tgnum_rhoimp = coefficient("tgnum_rhoimp")
heat_source_i = coefficient("heat_source_i")
ZK_i_perp_num_psin = coefficient("ZK_i_perp_num_psin")
tgnum_Ti = coefficient("tgnum_Ti")
visco_par_heating = coefficient("visco_par_heating")
implicit_heat_source = coefficient("implicit_heat_source")
Tie_min_neg = coefficient("Tie_min_neg")
aux_E0_Ti = coefficient("aux_E0_Ti")
heat_source_e = coefficient("heat_source_e")
ZK_e_perp_num_psin = coefficient("ZK_e_perp_num_psin")
tgnum_Te = coefficient("tgnum_Te")
ksi_ion_norm = coefficient("ksi_ion_norm")
aux_jre = coefficient("aux_jre")
aux_E0_Te = coefficient("aux_E0_Te")
power_dens_teleport_ju = coefficient("power_dens_teleport_ju")
heat_source_total = coefficient("heat_source")
ZK_perp_num_psin = coefficient("ZK_perp_num_psin")
tgnum_T = coefficient("tgnum_T")
T_min_neg = coefficient("T_min_neg")
aux_E0 = coefficient("aux_E0")

# -------------------------------------------------------------------------
# State-dependent work values
# -------------------------------------------------------------------------
# Quantities the element routine supplies together with the derivatives
# the Jacobian needs; declaring the interface lets the generator
# differentiate them by name.
eta = external_function(
    "eta600",
    arguments=("T", "rho", "rhoimp"),
    derivatives={"T": "deta_dT", "rho": "deta_dr0", "rhoimp": "deta_drimp0",},
    policy="piecewise_active",
    fortran_name="eta_T",
)

eta_num_T = external_function(
    "eta_num600",
    arguments=("T",),
    derivatives={"T": "deta_num_dT"},
    policy="piecewise_active",
    fortran_name="eta_num_T",
)

visco = external_function(
    "visco600",
    arguments=("T",),
    derivatives={"T": "dvisco_dT"},
    policy="piecewise_active",
    fortran_name="visco_T",
)

visco_num = external_function(
    "visco_num600",
    arguments=("T",),
    # The element routine computes dvisco_num_dT but deliberately does not
    # place that derivative in the momentum AMAT temperature columns.
    derivatives={},
    policy="frozen",
    fortran_name="visco_num_T",
)

alpha_i_state = external_function(
    "alpha_i600", arguments=(), derivatives={}, policy="frozen",
    fortran_name="alpha_i",
)

dvisco_state = external_function(
    "dvisco_state600", arguments=("T",),
    derivatives={"T": "d2visco_dT2"},
    policy="piecewise_active", fortran_name="dvisco_dT",
)

Sion_rate = external_function(
    "Sion_rate600", arguments=("T",), derivatives={"T": "dSion_dT"},
    policy="piecewise_active", fortran_name="Sion_T",
)

Srec_rate = external_function(
    "Srec_rate600", arguments=("T",), derivatives={"T": "dSrec_dT"},
    policy="piecewise_active", fortran_name="Srec_T",
)

alpha_e_state = external_function(
    "alpha_e600", arguments=("T",), derivatives={"T": "dalpha_e_dT"},
    policy="piecewise_active", fortran_name="alpha_e",
)

alpha_e_bis_state = external_function(
    "alpha_e_bis", arguments=("T",), derivatives={"T": "alpha_e_tri"},
    policy="piecewise_active", fortran_name="alpha_e_bis",
)

alpha_e_temperature = external_function(
    "alpha_e_temperature600", arguments=("T",), derivatives={"T": "alpha_e_bis"},
    policy="piecewise_active", fortran_name="alpha_e_T",
)

alpha_imp_bis_state = external_function(
    "alpha_imp_bis", arguments=("T",), derivatives={"T": "alpha_imp_tri"},
    policy="piecewise_active", fortran_name="alpha_imp_bis",
)

alpha_imp_temperature = external_function(
    "alpha_imp_temperature600", arguments=("T",),
    derivatives={"T": "alpha_imp_bis"},
    policy="piecewise_active", fortran_name="alpha_imp_T",
)

corr_neg_dens = external_function(
    "corr_neg_dens",
    arguments=("rho",),
    derivatives={"rho": "dr0_corr_dn"},
    policy="piecewise_active",
    fortran_name="corr_neg_dens",
)

W_dia_single = external_function(
    "W_dia_single600",
    arguments=("rho", "T"),
    derivatives={"rho": "W_dia_rho", "T": "W_dia_T"},
    policy="piecewise_active",
    fortran_name="W_dia",
)

W_dia_two = external_function(
    "W_dia_two600",
    arguments=("rho", "Ti"),
    derivatives={"rho": "W_dia_rho", "Ti": "W_dia_Ti",},
    policy="piecewise_active",
    fortran_name="W_dia",
)

ZKi_par = external_function(
    "ZKi_par600", arguments=("Ti",), derivatives={"Ti": "dZKi_par_dT"},
    policy="piecewise_active", fortran_name="ZKi_par_T",
)

ZKi_perp = external_function(
    "ZKi_prof600", arguments=("rho",), derivatives={"rho": "dZKi_prof_drho"},
    policy="piecewise_active", fortran_name="ZKi_prof",
)

visco_heating = external_function(
    "visco_heating600", arguments=("Te",),
    derivatives={"Te": "dvisco_dT_heating"},
    policy="piecewise_active", fortran_name="visco_T_heating",
)

Ti_e_exchange = external_function(
    "dTi_e600", arguments=("Ti", "Te", "rho", "rhoimp"),
    derivatives={"Ti": "ddTi_e_dTi", "Te": "ddTi_e_dTe", "rho": "ddTi_e_drho", "rhoimp": "ddTi_e_drhoimp"},
    policy="piecewise_active", fortran_name="dTi_e",
)

Ti_floor = external_function(
    "Ti_floor600", arguments=("Ti",), derivatives={"Ti": "dTi_floor"},
    policy="piecewise_active", fortran_name="Ti0_floor",
)

Ti_floor_exp = external_function(
    "Ti_floor_exp600", arguments=("Ti",),
    derivatives={"Ti": "dTi_floor_exp"},
    policy="piecewise_active", fortran_name="Ti_floor_exp",
)

corr_neg_dens_imp = external_function(
    "corr_neg_dens_imp", arguments=("rhoimp",),
    derivatives={"rhoimp": "drimp0_corr_dn"},
    policy="piecewise_active", fortran_name="corr_neg_dens_imp",
)

E_ion_bg_state = external_function(
    "E_ion_bg600", arguments=(), derivatives={}, policy="frozen",
    fortran_name="E_ion_bg",
)

ZKe_par = external_function(
    "ZKe_par600", arguments=("Te",), derivatives={"Te": "dZKe_par_dT"},
    policy="piecewise_active", fortran_name="ZKe_par_T",
)

ZKe_perp = external_function(
    "ZKe_prof600", arguments=("rho",), derivatives={"rho": "dZKe_prof_drho"},
    policy="piecewise_active", fortran_name="ZKe_prof",
)

eta_ohm_e = external_function(
    "eta_ohm600", arguments=("Te", "rho", "rhoimp"),
    derivatives={"Te": "deta_dT_ohm", "rho": "deta_dr0_ohm", "rhoimp": "deta_drimp0_ohm"},
    policy="piecewise_active", fortran_name="eta_T_ohm",
)

dE_ion_dT_state = external_function(
    "dE_ion_dT", arguments=("Te",), derivatives={}, policy="frozen",
)

E_ion_state = external_function(
    "E_ion600", arguments=("Te",), derivatives={"Te": "dE_ion_dT"},
    policy="piecewise_active", fortran_name="E_ion",
)

LradDrays = external_function(
    "LradDrays600", arguments=("Te",), derivatives={"Te": "dLradDrays_dT"},
    policy="piecewise_active", fortran_name="LradDrays_T",
)

LradDcont = external_function(
    "LradDcont600", arguments=("Te",),
    derivatives={"Te": "dLradDcont_dT_corr"},
    policy="piecewise_active", fortran_name="LradDcont_corr",
)

frad_bg_state = external_function(
    "frad_bg600", arguments=("Te",), derivatives={"Te": "dfrad_bg_dT"},
    policy="piecewise_active", fortran_name="frad_bg",
)

Lrad_state = external_function(
    "Lrad600", arguments=("Te",), derivatives={"Te": "dLrad_dT"},
    policy="piecewise_active", fortran_name="Lrad",
)

Te_i_exchange = external_function(
    "dTe_i600", arguments=("Ti", "Te", "rho", "rhoimp"),
    derivatives={"Ti": "ddTe_i_dTi", "Te": "ddTe_i_dTe", "rho": "ddTe_i_drho", "rhoimp": "ddTe_i_drhoimp"},
    policy="piecewise_active", fortran_name="dTe_i",
)

Te_floor = external_function(
    "Te_floor600", arguments=("Te",), derivatives={"Te": "dTe_floor"},
    policy="piecewise_active", fortran_name="Te0_floor",
)

Te_floor_exp = external_function(
    "Te_floor_exp600", arguments=("Te",),
    derivatives={"Te": "dTe_floor_exp"},
    policy="piecewise_active", fortran_name="Te_floor_exp",
)

corr_neg_dens_n = external_function(
    "corr_neg_dens_n", arguments=("rhon",),
    derivatives={"rhon": "drn0_corr_dn"},
    policy="piecewise_active", fortran_name="corr_neg_dens_n",
)

ZK_par = external_function(
    "ZK_par600", arguments=("T",), derivatives={"T": "dZK_par_dT"},
    policy="piecewise_active", fortran_name="ZK_par_T",
)

ZK_perp = external_function(
    "ZK_prof600", arguments=("rho",), derivatives={"rho": "dZK_prof_drho"},
    policy="piecewise_active", fortran_name="ZK_prof",
)

T_floor = external_function(
    "T_floor600", arguments=("T",), derivatives={"T": "dT_floor"},
    policy="piecewise_active", fortran_name="T0_floor",
)

T_floor_exp = external_function(
    "T_floor_exp600", arguments=("T",), derivatives={"T": "dT_floor_exp"},
    policy="piecewise_active", fortran_name="T_floor_exp",
)

# -------------------------------------------------------------------------
# Operators and shared term groups
# -------------------------------------------------------------------------
def _rho_hat(value=rho):
    return R**2 * value

def _B2(flux):
    """Total magnetic field squared (B^2)"""
    return (F0**2 + grad(flux)[0]**2 + grad(flux)[1]**2) / R**2

def _B_dot_grad(value, flux=psi, *, st_form=False):
    """Return B.grad(value). st_form computes the Poiss bracket in st coordinates"""

    if st_form:
        return (F0 / R * dphi(value) + poiss_bracket_st(value, flux) / xjac) / R

    return (F0 / R * dphi(value) + poiss_bracket(value, flux)) / R

def _diffusion_tot_intg_by_parts(test, value):
    """Total diffusion operator, including its toroidal part, integrated by parts."""

    return dot(grad(test), grad(value)) + dphi(test) * dphi(value) / R**2

def _u_convection(density, *, bracket=poiss_bracket_st):
    """Advection and compression of a particle density by u flow."""

    return R * bracket(density, u) / xjac + 2 * density * dZ(u) 

def _parallel_convection(quantity, vpar=vpar):
    """Convection of a quantity density along the field by ``vpar``."""

    return (
        - vpar     * _B_dot_grad(quantity, st_form=True)
        - quantity * _B_dot_grad(vpar,    st_form=True) 
    )

def _dens_tgnum_intg_by_parts(v, density):
    """Taylor-Galerkin stabilization of a particle-density equation."""

    return (
        -sp.Rational(1, 4) * R**2 * poiss_bracket(density, u) * poiss_bracket(v, u) * tstep
        -sp.Rational(1, 4) * vpar**2 * _B_dot_grad(density) * _B_dot_grad(v) * tstep
    )

def _par_diff_intg_by_parts(v, density):
    """Parallel diffusion term along the field integrated by parts."""

    return - _B_dot_grad(v) * _B_dot_grad(density) / _B2(psi) 

def _diamagnetic_pressure(*, with_TiTe=False):
    """Ion pressure ``Pi = (rho + alpha_i*rhoimp)*Ti`` of ``construct_pressure``.

    The one-temperature model evolves the total temperature and supplies
    Ti0 = T0/2, so the ion temperature is half of the evolved field there.
    ``alpha_i`` carries no temperature dependence in either model.
    """

    ion_temperature = Ti if with_TiTe else T / 2
    return (rho + rhoimp * alpha_i_state()) * ion_temperature

def _diamagnetic_viscosity(*, with_TiTe=False):
    """The element routine's ``W_dia``.

    It is built from the ion pressure alone, so it carries no explicit
    electron-temperature dependence.
    """

    return W_dia_two(rho, Ti) if with_TiTe else W_dia_single(rho, T)


def _press_u_convection(pressure, *, bracket):
    """Convection and compression of a pressure by the u (ExB) flow."""

    return (R * bracket(pressure, u) / xjac + 2 * GAMMA * pressure * dZ(u))


def _press_parallel_convection(pressure):
    """Convection and compression of a pressure along the field by ``vpar``."""

    return (
        - vpar * _B_dot_grad(pressure, st_form=True)
        - GAMMA * pressure * _B_dot_grad(vpar, st_form=True)
    )


def _energy_tgnum(v, pressure, factor):
    """Taylor-Galerkin stabilization of an energy equation.   `."""

    return (
        -factor * sp.Rational(1, 4) * R**2 * poiss_bracket(pressure, u) * poiss_bracket(v, u) * tstep
        - factor * sp.Rational(1, 4) * vpar**2 * _B_dot_grad(pressure) * _B_dot_grad(v) * tstep
    )


def _heating_floor(v, exponential, floor, minimum):
    """Implicit heating that keeps a temperature away from its floor.   `."""

    return implicit_heat_source * (gamma - 1) * v * (sp.Rational(1,2) * minimum * (1 + exponential) - floor)


def _released_kinetic_energy(T_or_Te):
    """Particle sources whose kinetic energy is released into the ions."""

    return (
        (rho + alpha_e_state(T_or_Te) * rhoimp) * rhon * Sion_rate(T_or_Te)
        + particle_source + source_pellet + source_bg_drift + source_imp_drift
    )


def _parallel_viscous_heating(v):
    """Heating from the parallel (``vpar``) viscosity. """

    grad_vpar = grad(vpar)
    return (GAMMA - 1) * visco_par_heating * (v * dot(grad_vpar, grad_vpar) + vpar * dot(grad(v), grad_vpar))


def _u_viscous_heating(v, heating):
    """Heating from the perpendicular (poloidal ``u`` flow) viscosity."""

    return (
        - (GAMMA - 1) * v * heating * R**2 * visco_fact_old  * dot(grad(u), grad(omega))
        - (GAMMA - 1) * v * heating * 2 * R * visco_fact_new * omega * dR(u)
        - (GAMMA - 1) * v * heating * visco_fact_new * (dR(u) * dR(dphi(dphi(u))) + dZ(u) * dZ(dphi(dphi(u))))
    )


def _kinetic_coupling(v):
    """Energy and momentum handed over by the kinetic neutral/impurity model"""

    return (
        (gamma - 1) * sp.Rational(1, 2) * v * aux_rho0 * vpar**2 * _B2(psi)
        - (gamma - 1) * v * aux_mom_par0 * vpar
    )


def _ionization_energy_transport(v, energy, temperature):
    """Transport and diffusive flux of the ionization potential energy."""

    d_par_tot = D_par_local + D_par_sc_num * tau_sc - D_prof
    d_par_imp_tot = D_par_local_imp + D_par_imp_sc_num * tau_sc - D_prof_imp
    bb2 = _B2(psi)
    return (GAMMA - 1) * (
        + v * R * poiss_bracket_st(energy, u) / xjac
        + 2 * v * energy * dZ(u)
        - v * F0 * vpar * dphi(energy) / R**2
        - v * vpar * poiss_bracket_st(energy, psi) / (R * xjac)
        - v * energy * _B_dot_grad(vpar, st_form=True)
        - E_ion_state(temperature) * d_par_imp_tot / bb2 * _B_dot_grad(v) * _B_dot_grad(rhoimp)
        - E_ion_state(temperature) * D_prof_imp * _diffusion_tot_intg_by_parts(v, rhoimp)
        - E_ion_bg_state() * d_par_tot / bb2 * _B_dot_grad(v) * _B_dot_grad(rho - rhoimp)
        - E_ion_bg_state() * D_prof * _diffusion_tot_intg_by_parts(v, rho - rhoimp)
    )


# -------------------------------------------------------------------------
# Equations
# -------------------------------------------------------------------------
# In the order of the element routine.
def induction_equation_1( with_TiTe: bool = False):
    """Weak form of the induction equation for var_psi."""

    v = test_function("v")
    T_or_Te = Te if with_TiTe else T
    Te_gen = Te if with_TiTe else T / 2
    # Electron pressure
    Pe = rho * Te_gen + rhoimp * alpha_e_temperature(Te_gen)
    rho_corr = corr_neg_dens(rho)
    dV = R * xjac

    # Time derivative
    A = v * psi / R**2 * dV   

    # RHS
    B = (
        # -B.grad u
        + v * poiss_bracket_st(psi, u) / (R * xjac) - v * F0 / R**2 * dphi(u)

        # eta*j
        + v * eta(T_or_Te, rho, rhoimp)*(j - coefficient("current_source") - Jb) / R**2

        # hyper-resistivity
        + eta_num_T(T_or_Te) * dot(grad(v), grad(j)) / R

        # Runaway current coupling, -eta*j_RE
        -v * eta(T_or_Te, rho, rhoimp) * aux_jre_ind / R**2

        # B.grad Pe diamagnetic term (no impurity contribution)
        - v * 2 * tauIC / (rho_corr * _B2(psi)) * F0**2 / R**3 * poiss_bracket_st(psi, Pe) / xjac
        + v * 2 * tauIC / (rho_corr * _B2(psi)) * F0**3 / R**4 * dphi(Pe) 
        
    ) * dV

    return EvolutionEquation("model600_induction", v, A, B)

def momentum_equation_2(*, with_TiTe=False, include_neo=False):
    """Model-600 perpendicular momentum equation (``var_u``). """

    v = test_function("v")
    T_or_Te = Te if with_TiTe else T
    alpha_e_value = alpha_e_state(T_or_Te)
    sion_rate = Sion_rate(T_or_Te)
    srec_rate = Srec_rate(T_or_Te)
    dV = R * xjac

    # In the two-temperature branch p0 is the sum of ion and electron
    # pressures; the single-temperature branch uses the unified T field.
    if with_TiTe:
        pressure = rho * (Ti + Te) + rhoimp * (alpha_i_state() * Ti + alpha_e_temperature(Te))
    else:
        pressure = rho * T + rhoimp * alpha_imp_temperature(T)

    rho_hat = _rho_hat()
    grad_v = grad(v)
    grad_u = grad(u)
    grad_omega = grad(omega)
    lap_v = laplacian(v)
    lap_omega = laplacian(omega)
    pi = _diamagnetic_pressure(with_TiTe=with_TiTe)
    W_dia = _diamagnetic_viscosity(with_TiTe=with_TiTe)
    ne = rho + alpha_e_value * rhoimp
    neutral_sources = ne * rhon * sion_rate - ne * (rho - rhoimp) * srec_rate

    # Time derivative
    A = R**2 * dV * (
        - corr_neg_dens(freeze(rho)) * dot(grad_v, grad_u)
        - fact_conservative_u * rho * dot(grad_v, grad(freeze(u))) # Not sure one should freeze u
    )

    # The RHS, in the order the element routine builds it.
    B = (
        # Perpendicular inertia (poloidal Jacobian of the kinetic energy).
        - R * sp.Rational(1, 2) * (dR(u)**2 + dZ(u)**2) * poiss_bracket(v, rho_hat)

        # Vorticity advection.
        - R**3 * rho * omega * poiss_bracket_st(v, u) / xjac

        # jxB force term. which here is B.gra j
        + v * poiss_bracket_st(psi, j) / dV
        - v * F0 * dphi(j) / R**2

        # Pressure advection.
        + R * poiss_bracket_st(v, pressure) / xjac

        # Perpendicular and toroidal viscosity.
        - visco(T_or_Te) * R**2 * visco_fact_old * dot(grad_v, grad_omega)
        - 2 * visco(T_or_Te) * R * visco_fact_new * omega * dR(v)
        - visco(T_or_Te) * visco_fact_new * (dR(v) * dR(dphi(dphi(u))) + dZ(v) * dZ(dphi(dphi(u))))
        - visco_num(T_or_Te) * lap_v * lap_omega / R

        # Diamagnetic pressure advection and diamagnetic viscosity.
        - v * tauIC * 2 * R**3 * poiss_bracket_st(pi, omega) / xjac
        - tauIC * 2 * R**2 * dZ(pi) * dot(grad_v, grad_u)
        - v * tauIC * 2 * R**3 * (dR(dZ(u)) * (dR(dR(pi)) - dZ(dZ(pi))) - dR(dZ(pi)) * (dR(dR(u)) - dZ(dZ(u))))
        + dvisco_state(T_or_Te) * W_dia * dot(grad(Ti if with_TiTe else T / 2), grad_v)
        + visco(T_or_Te) * W_dia * lap_v

        # Conservative form of the momentum equation.
        + fact_conservative_u * dot(grad_v, grad_u) * (
            - R * poiss_bracket(rho_hat, u)
            + F0 * (rho * dphi(vpar) + vpar * dphi(rho))
            + R * rho * poiss_bracket(vpar, psi)
            + R * vpar * poiss_bracket(rho, psi)
        )

        # Runaway/auxiliary pressure contributions
        - v * (aux_P_par_re + aux_P_perp_re)
        + R * (-aux_divPIR_perp * dZ(v) + aux_divPIZ_perp * dR(v))

        # Momentum carried by the neutral and particle sources
        + R**2 * dot(grad_v, grad_u) * (
            (1 - delta_n_convection) * neutral_sources
            + (1 - fact_conservative_u) * (particle_source + source_pellet + source_bg_drift + source_imp_drift)
        )

        # Taylor-Galerkin stabilization. 
        - tgnum_u * sp.Rational(1, 4) * R**4 * rho* poiss_bracket(omega, u) * poiss_bracket(v, u) * tstep
        - tgnum_u * sp.Rational(1, 4) * R**2 * omega * fact_conservative_u * poiss_bracket(rho_hat, u) * poiss_bracket(v, u) * tstep

    ) * dV

    # Neo-classical terms still missing!

    return EvolutionEquation("model600_momentum", v, A, B, kind="evolution", amat_variation_overrides={})

def current_constraint_equation_zj():
    """Model-600 current-definition equation for ``var_zj``."""

    v = test_function("v")
    C = (dot(grad(v), grad(psi)) + v * j) / R * xjac
    return ConstraintEquation("model600_current_constraint", v, C, kind="static")

def vorticity_constraint_equation_w():
    """Model-600 vorticity-definition equation for ``var_w``."""

    v = test_function("v")
    C = (dot(grad(v), grad(u)) + v * omega) * R * xjac
    return ConstraintEquation("model600_vorticity_constraint", v, C, kind="static")

def density_equation_rho(*, with_TiTe=False):
    """Model-600 density equation (``var_rho``)."""

    v = test_function("v")
    T_or_Te = Te if with_TiTe else T
    rho_main = rho - rhoimp
    ne = rho + alpha_e_state(T_or_Te) * rhoimp
    D_par_tot = D_par_local + D_par_sc_num * tau_sc 
    D_par_imp_tot = D_par_local_imp + D_par_imp_sc_num * tau_sc 
    psi_gradient = grad(psi)
    dV = R * xjac

    # Time derivative
    A = v * rho * dV

    # Other terms and RHS
    B = (
        # particle sources
        + v * (particle_source+source_pellet+source_bg_drift+source_imp_drift+aux_rho0) 

        # sources/sinks due to neutrals
        + v * ne * rhon * Sion_rate(T_or_Te) 
        - v * ne * rho_main * Srec_rate(T_or_Te) 

        # convection/compression by u variable (ExB)
        + v * _u_convection(rho) 

        # parallel convection/compression by vpar
        + v * _parallel_convection(rho)  

        # parallel and perpendicular diffusion of main ions
        + (D_par_tot-D_prof) * _par_diff_intg_by_parts(v, rho_main) 
        - D_prof        * _diffusion_tot_intg_by_parts(v, rho_main) 

        # parallel and perpendicular diffusion of impurities (total mass)
        + (D_par_imp_tot-D_prof_imp) * _par_diff_intg_by_parts(v, rhoimp)   
        - D_prof_imp *     _diffusion_tot_intg_by_parts(v, rhoimp)

        # numerical stabilization
        + tgnum_rho * _dens_tgnum_intg_by_parts(v, rho) 
        - D_perp_num_psin * laplacian(v) * laplacian(rho) 

        # Diamagnetic drift
        + v * 2 * tauIC * 2 * dZ(_diamagnetic_pressure(with_TiTe=with_TiTe))

        # Pinch term, not sure about this one here
        - V_prof_pinch / sp.sqrt(psi_gradient[0]**2 + psi_gradient[1]**2) * dot(grad(v), psi_gradient) * rho 

    ) * dV
    return EvolutionEquation("model600_density", v, A, B)

def parallel_velocity_equation_vpar(*, with_TiTe=False, st_form=True):
    """Model-600 parallel velocity equation (``var_vpar``)."""

    v = test_function("v")
    T_or_Te = Te if with_TiTe else T
    bb2 = _B2(psi)
    dV = R * xjac
    rho_hat = _rho_hat()
    psi_gradient = grad(psi)
    ne = rho + alpha_e_state(T_or_Te) * rhoimp
    pressure = rho * (Ti + Te) if with_TiTe else rho * T
    if with_TiTe:
        pressure += rhoimp * (alpha_i_state() * Ti + alpha_e_temperature(Te))
    else:
        pressure += rhoimp * alpha_imp_temperature(T)

    visco_par_eff = visco_par + visco_par_sc_num * tau_sc
    source_dens_tot = particle_source + source_pellet + source_bg_drift + source_imp_drift

    Bdot_grad = functools.partial(_B_dot_grad, st_form=st_form)
    # Only the flux is varied in the prescribed rotation profile, so
    # ``grad(Vt) = dV_dpsi_source*grad(psi)`` is the faithful form.
    rotation_shear = tuple(
        parallel - dV_dpsi_source * flux
        for parallel, flux in zip(grad(vpar), grad(psi))
    )

    # Time derivative
    # The parallel momentum density is rho*v_par*B**2. The element routinefreezes the density correction in the tangent..
    A = v * (
        corr_neg_dens(freeze(rho)) * vpar * bb2
        + fact_conservative_u * rho * freeze(vpar) * bb2
    ) * dV

    # The RHS and others
    B = (
        # Parallel pressure gradient.
        - v * _B_dot_grad(pressure, st_form=True) 

        # Parallel advection of the kinetic energy, 0.5*v_par**2*B**2.
        + sp.Rational(1, 2) * vpar**2 * bb2 * (rho * Bdot_grad(v) + v * Bdot_grad(rho))

        # Term to obtain conservative form of momentum equation
        + fact_conservative_u * v * vpar * bb2 * (poiss_bracket(rho_hat, u)/R - _B_dot_grad(rho * vpar))

        # Numerical and physical parallel viscosities.
        - visco_par_num * laplacian(v) * laplacian(vpar) 
        - visco_par_par * F0**2 / (R**2 * bb2) * _B_dot_grad(vpar) * _B_dot_grad(v)
        - visco_par_eff * dot(grad(v), rotation_shear)

        # External momentum sources
        + v * aux_mom_par0 

        # Momentum carried by the particle sources; not active in conservative form
        - v * source_dens_tot * vpar * bb2 * (1 - fact_conservative_u) 
        - v * aux_rho0       *  vpar * bb2 * (1 - fact_conservative_u)

        # This one should probably be multiplied by (1 - fact_conservative_u
        + (1 - delta_n_convection) * v * ne * vpar * bb2 * ((rho - rhoimp) * Srec_rate(T_or_Te) -  rhon * Sion_rate(T_or_Te))

        # Taylor-Galerkin stabilization.
        - tgnum_vpar * sp.Rational(1, 4) * rho * vpar**2 * bb2* _B_dot_grad(vpar, st_form=True) * _B_dot_grad(v, st_form=True) * tstep
        - tgnum_vpar * sp.Rational(1, 4) * v * vpar**2 * bb2 * (1 - fact_conservative_u) * _B_dot_grad(vpar, st_form=True) * _B_dot_grad(rho, st_form=True) * tstep
        - tgnum_vpar * sp.Rational(1, 4) * vpar**3 * bb2 * fact_conservative_u * _B_dot_grad(rho, st_form=True) * _B_dot_grad(v, st_form=True) * tstep

        # Not sure this pinch term should be here
        + V_prof_pinch / sp.sqrt(psi_gradient[0]**2 + psi_gradient[1]**2) * dot(psi_gradient, grad(vpar)) * rho * v

    ) * dV

    return EvolutionEquation("model600_parallel_velocity", v, A, B)

def impurity_density_equation_rhoimp(*, with_TiTe=False):
    """Model-600 impurity-density equation (``var_rhoimp``)"""

    v = test_function("v")
    D_par_imp_tot = D_par_local_imp + D_par_imp_sc_num * tau_sc
    dV = R * xjac

    # Time derivative
    A = v * rhoimp * dV

    # Other terms and RHS
    B = (
        # particle sources
        + v * source_imp_drift 

        # convection/compression by u variable (ExB)
        + v * _u_convection(rhoimp) 

        # parallel convection/compression by vpar
        + v * _parallel_convection(rhoimp)

        # parallel and perpendicular diffusion
        + (D_par_imp_tot-D_prof_imp) * _par_diff_intg_by_parts(v, rhoimp)
        - D_prof_imp  * _diffusion_tot_intg_by_parts(v, rhoimp)

        # numerical stabilization
        + tgnum_rhoimp * _dens_tgnum_intg_by_parts(v, rhoimp)
        - Dn_perp_num * laplacian(v) * laplacian(rhoimp)
    ) * dV
    
    return EvolutionEquation("model600_impurity_density", v, A, B)

def ion_energy_equation_Ti(*, st_form=True):
    """Model-600 ion energy equation (``var_Ti``).

    ``st_form`` selects the spelling of the poloidal advection bracket.
    """

    v = test_function("v")
    bracket = poiss_bracket_st if st_form else (lambda a, b: poiss_bracket(a, b) * xjac)
    dV = R * xjac
    ion_density = rho + alpha_i_state() * rhoimp
    ion_pressure = ion_density * Ti

    # Time derivative
    A = v * (
        corr_neg_dens(rho) + alpha_i_state() * corr_neg_dens_imp(rhoimp)
    ) * Ti * dV

    # RHS, many terms are common to T/Ti/Te (so many operators are reused to avoid duplications)
    B = (
        # convection/compression by u (ExB)
        + v * _press_u_convection(ion_pressure, bracket=bracket)

        # parallel convection/compression by vpar
        + v * _press_parallel_convection(ion_pressure)

        # parallel heat conduction
        + (ZKi_par(Ti) - ZKi_perp(rho)) * _par_diff_intg_by_parts(v, Ti)

        # perpendicular (total) heat conduction
        - ZKi_perp(rho) * _diffusion_tot_intg_by_parts(v, Ti)

        # heat sources/sinks
        + v * (heat_source_i + aux_E0_Ti)

        # energy exchange with electrons
        + v * Ti_e_exchange(Ti, Te, rho, rhoimp)

        # neutral recombination sink
        - v * Ti * corr_neg_dens(rho)**2 * Srec_rate(Te)

        # kinetic energy released by particle sources (friction heating)
        + v * (GAMMA - 1) / 2 * (vpar**2 * _B2(psi) + R**2 * (dR(u)**2 + dZ(u)**2)) * _released_kinetic_energy(Te)

        # parallel viscous heating
        + _parallel_viscous_heating(v)

        # perpendicular (u) viscous heating
        + _u_viscous_heating(v, visco_heating(Te))

        # energy and momentum from the kinetic neutral/impurity model
        + _kinetic_coupling(v)

        # numerical stabilization
        + _energy_tgnum(v, ion_pressure, tgnum_Ti)
        - ZK_i_perp_num_psin * laplacian(v) * laplacian(Ti)

        # implicit heat source to avoid negative temperatures
        + _heating_floor(v, Ti_floor_exp(Ti), Ti_floor(Ti), Tie_min_neg)

    ) * dV
    return EvolutionEquation("model600_ion_energy", v, A, B)

def electron_energy_equation_Te(*, st_form=True):
    """Model-600 electron energy equation (``var_Te``).

    ``st_form`` applies to the electron-pressure advection only; the
    ionization-energy advection is always written in element coordinates.
    """

    v = test_function("v")
    bracket = poiss_bracket_st if st_form else (lambda a, b: poiss_bracket(a, b) * xjac)
    dV = R * xjac
    electron_pressure = rho * Te + rhoimp * alpha_e_temperature(Te)
    electron_density = corr_neg_dens(rho) + alpha_e_state(Te) * corr_neg_dens_imp(rhoimp)
    ionization_energy = E_ion_state(Te) * rhoimp + E_ion_bg_state() * (rho - rhoimp)

    # Time derivative
    A = v * (
        corr_neg_dens(rho) * Te + corr_neg_dens_imp(rhoimp) * alpha_e_temperature(Te)
        + (GAMMA - 1) * ionization_energy
    ) * dV

    # RHS, many terms are common to T/Ti/Te (so many operators are reused to avoid duplications)
    B = (
        # convection/compression by u (ExB)
        + v * _press_u_convection(electron_pressure, bracket=bracket)

        # parallel convection/compression by vpar
        + v * _press_parallel_convection(electron_pressure)

        # parallel heat conduction
        + (ZKe_par(Te) - ZKe_perp(rho)) * _par_diff_intg_by_parts(v, Te)

        # perpendicular (total) heat conduction
        - ZKe_perp(rho) * _diffusion_tot_intg_by_parts(v, Te)

        # heat sources/sinks
        + v * ( heat_source_e + aux_E0_Te + power_dens_teleport_ju)

        # Ionization sink due to neutrals
        - v * ksi_ion_norm * (rho + alpha_e_state(Te) * rhoimp) * rhon * Sion_rate(Te)

        # energy exchange with ions
        + v * Te_i_exchange(Ti, Te, rho, rhoimp)

        # ohmic heating
        + v * (GAMMA - 1) * eta_ohm_e(Te, rho, rhoimp) * ((j - aux_jre) / R)**2

        # neutral (line) radiation
        - v * electron_density * corr_neg_dens_n(rhon) * LradDrays(Te)

        # background radiation
        - v * electron_density * (corr_neg_dens(rho) - corr_neg_dens_imp(rhoimp)) * LradDcont(Te)
        - v * electron_density * frad_bg_state(Te)

        # impurity radiation
        - v * electron_density * corr_neg_dens_imp(rhoimp) * Lrad_state(Te)

        # ionization potential transport
        + _ionization_energy_transport(v, ionization_energy, Te)

        # numerical stabilization
        + _energy_tgnum(v, electron_pressure, tgnum_Te)
        - ZK_e_perp_num_psin * laplacian(v) * laplacian(Te)
        
        # implicit heat source to avoid negative temperatures
        + _heating_floor(v, Te_floor_exp(Te), Te_floor(Te), Tie_min_neg)

    ) * dV
    return EvolutionEquation("model600_electron_energy", v, A, B)

def total_energy_equation_T(*, st_form=True):
    """Model-600 single-temperature energy equation (``var_T``)"""

    v = test_function("v")
    bracket = poiss_bracket_st if st_form else (lambda a, b: poiss_bracket(a, b) * xjac)
    dV = R * xjac
    pressure = rho * T + rhoimp * alpha_imp_temperature(T)
    electron_density = corr_neg_dens(rho) + alpha_e_state(T) * corr_neg_dens_imp(rhoimp)
    ionization_energy = E_ion_state(T) * rhoimp + E_ion_bg_state() * (rho - rhoimp)

    # Time derivative
    A = v * (
        corr_neg_dens(rho) * T + corr_neg_dens_imp(rhoimp) * alpha_imp_temperature(T)
        + (GAMMA - 1) * ionization_energy
    ) * dV

    # RHS, many terms are common to T/Ti/Te (so many operators are reused to avoid duplications)
    B = (
        # convection/compression by u (ExB)
        + v * _press_u_convection(pressure, bracket=bracket)

        # parallel convection/compression by vpar
        + v * _press_parallel_convection(pressure)

        # parallel heat conduction
        + (ZK_par(T) - ZK_perp(rho)) * _par_diff_intg_by_parts(v, T)

        # perpendicular (total) heat conduction
        - ZK_perp(rho) * _diffusion_tot_intg_by_parts(v, T)

        # heat sources/sinks
        + v * (heat_source_total + aux_E0 + power_dens_teleport_ju)

        # Ionization sink due to neutrals
        - v * ksi_ion_norm * (rho + alpha_e_state(T) * rhoimp) * rhon * Sion_rate(T)

        # recombination sink due to neutrals
        - (GAMMA - 1) * v * sp.Rational(1, 2) * T * corr_neg_dens(rho)**2 * Srec_rate(T)

        # ohmic heating
        + v * (GAMMA - 1) * eta_ohm_e(T, rho, rhoimp) * ((j - aux_jre) / R)**2

        # neutral (line) radiation
        - v * electron_density * corr_neg_dens_n(rhon) * LradDrays(T)

        # background radiation
        - v * electron_density * (corr_neg_dens(rho) - corr_neg_dens_imp(rhoimp)) * LradDcont(T)
        - v * electron_density * frad_bg_state(T)

        # impurity radiation
        - v * electron_density * corr_neg_dens_imp(rhoimp) * Lrad_state(T)

        # ionization potential transport
        + _ionization_energy_transport(v, ionization_energy, T)

        # kinetic energy released by particle sources (friction heating)
        + v * (GAMMA - 1) / 2 * (vpar**2 * _B2(psi) + R**2 * (dR(u)**2 + dZ(u)**2)) * _released_kinetic_energy(T)

        # parallel viscous heating
        + _parallel_viscous_heating(v)

        # perpendicular (u) viscous heating
        + _u_viscous_heating(v, visco_heating(T))

        # energy and momentum from the kinetic neutral/impurity model
        + _kinetic_coupling(v)

        # numerical stabilization
        + _energy_tgnum(v, pressure, tgnum_T)
        - ZK_perp_num_psin * laplacian(v) * laplacian(T)

        # implicit heat source to avoid negative temperatures
        + _heating_floor(v, T_floor_exp(T), T_floor(T), T_min_neg)

    ) * dV
    return EvolutionEquation("model600_total_energy", v, A, B)


