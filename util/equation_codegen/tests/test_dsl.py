"""The symbolic DSL: fields and roles, operators, externals, sign conventions."""

import sys
import unittest
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from jorek_equations import (  # noqa: E402
    ConstraintEquation,
    EvolutionEquation,
    FieldRole,
    MissingDerivativeError,
    bracket,
    coefficient,
    dR,
    dZ,
    dphi,
    poiss_bracket_st,
    expand_brackets,
    expand_derivatives,
    external_function,
    field,
    fortran,
    freeze,
    grad,
    dot,
    test_function,
    variation,
)
from jorek_equations.model199 import FIELDS, eta, j, psi, u



class SymbolsAndOperatorsTest(unittest.TestCase):
    def test_model199_field_order_and_fortran_names(self):
        self.assertEqual(
            [value.name for value in FIELDS],
            ["psi", "u", "j", "omega", "rho", "T"],
        )
        self.assertEqual(FIELDS[2].fortran_current, "zj0")
        self.assertEqual(FIELDS[3].fortran_trial, "w")

    def test_roles_are_distinct(self):
        psi = field("psi_role_test")
        self.assertEqual(psi.current.role, FieldRole.CURRENT)
        self.assertEqual(psi.increment.role, FieldRole.DELTA)
        self.assertEqual(psi.previous_increment.role, FieldRole.PREVIOUS_DELTA)
        self.assertEqual(psi.trial.role, FieldRole.TRIAL)
        self.assertNotEqual(psi.current, psi.trial)

    def test_bare_field_is_abstract_until_linearization(self):
        psi = field("psi_abstract_test")
        u = field("u_abstract_test")
        expression = psi * u

        self.assertEqual(variation(expression, psi), psi.increment * u.current)
        self.assertEqual(variation(expression, u), psi.current * u.increment)

    def test_derivative_variation_commutes(self):
        psi = field("psi_derivative_test")
        actual = variation(dR(psi.current), psi)
        self.assertEqual(actual, dR(psi.increment))

    def test_bracket_product_rule(self):
        psi = field("psi_bracket_test")
        u = field("u_bracket_test")
        expression = bracket(psi.current, u.current)

        self.assertEqual(
            variation(expression, psi, direction=FieldRole.TRIAL),
            bracket(psi.trial, u.current),
        )
        self.assertEqual(
            variation(expression, u, direction=FieldRole.TRIAL),
            bracket(psi.current, u.trial),
        )

    def test_gradient_dot_product_rule(self):
        u = field("u_gradient_test")
        expression = sum(component**2 for component in grad(u.current))
        expected = 2 * (
            dR(u.current) * dR(u.trial)
            + dZ(u.current) * dZ(u.trial)
        )
        actual = variation(expression, u, direction=FieldRole.TRIAL)
        self.assertEqual(sp.expand(actual - expected), 0)

    def test_frozen_expression_has_zero_variation(self):
        rho = field("rho_frozen_test")
        u = field("u_frozen_test")
        expression = freeze(rho.current) * dR(u.current)

        self.assertEqual(variation(expression, rho), 0)
        self.assertEqual(
            variation(expression, u),
            freeze(rho.current) * dR(u.increment),
        )

    def test_test_function_and_coefficients_do_not_vary(self):
        psi = field("psi_test_function_test")
        v = test_function("v_test_function_test")
        Rloc = coefficient("R_test_function_test")
        expression = v * psi.current / Rloc

        expected = v * psi.trial / Rloc
        actual = variation(expression, psi, direction=FieldRole.TRIAL)
        self.assertEqual(actual, expected)



class ExternalFunctionsTest(unittest.TestCase):
    def test_active_function_uses_supplied_derivative(self):
        T = field("T_external_active_test")
        j = field("j_external_active_test")
        eta = external_function(
            "eta_external_active_test",
            arguments=("T",),
            derivatives={"T": "deta_dT_external_active_test"},
        )
        expression = eta(T.current) * j.current

        expected_T = (
            sp.Function("deta_dT_external_active_test")(T.current)
            * T.trial
            * j.current
        )
        expected_j = eta(T.current) * j.trial

        self.assertEqual(
            variation(expression, T, direction=FieldRole.TRIAL), expected_T
        )
        self.assertEqual(
            variation(expression, j, direction=FieldRole.TRIAL), expected_j
        )

    def test_piecewise_active_uses_branch_derivative_interface(self):
        T = field("T_external_piecewise_test")
        eta = external_function(
            "eta_external_piecewise_test",
            arguments=("T",),
            derivatives={"T": "deta_dT_external_piecewise_test"},
            policy="piecewise_active",
        )

        expected = (
            sp.Function("deta_dT_external_piecewise_test")(T.current)
            * T.increment
        )
        self.assertEqual(variation(eta(T.current), T), expected)

    def test_frozen_external_function_does_not_propagate_dependency(self):
        psi = field("psi_external_frozen_test")
        Dperp = external_function(
            "Dperp_external_frozen_test",
            arguments=("psi_norm",),
            derivatives={},
            policy="frozen",
        )

        self.assertEqual(variation(Dperp(psi.current), psi), 0)

    def test_missing_active_derivative_is_an_error(self):
        psi = field("psi_external_missing_test")
        profile = external_function(
            "profile_external_missing_test",
            arguments=("psi",),
            derivatives={},
            policy="active",
        )

        with self.assertRaises(MissingDerivativeError):
            variation(profile(psi.current), psi)



class EquationConventionTest(unittest.TestCase):
    def test_evolution_rhs_and_amat_signs(self):
        psi = field("psi_evolution_test")
        j = field("j_evolution_test")
        T = field("T_evolution_test")
        v = test_function("v_evolution_test")
        Rloc = coefficient("R_evolution_test")
        timestep = coefficient("tstep_evolution_test")
        theta = coefficient("theta_evolution_test")
        zeta = coefficient("zeta_evolution_test")
        eta = external_function(
            "eta_evolution_test",
            arguments=("T",),
            derivatives={"T": "deta_dT_evolution_test"},
        )

        A = v * psi.current / Rloc
        B = v * eta(T.current) * j.current / Rloc
        result = EvolutionEquation("induction_test", v, A, B).linearize(
            fields=(psi, j, T),
            timestep=timestep,
            theta=theta,
            zeta=zeta,
        )

        expected_rhs = (
            timestep * B
            + zeta * v * psi.previous_increment / Rloc
        )
        expected_psi = (1 + zeta) * v * psi.trial / Rloc
        expected_j = -theta * timestep * v * eta(T.current) * j.trial / Rloc
        expected_T = (
            -theta
            * timestep
            * v
            * sp.Function("deta_dT_evolution_test")(T.current)
            * T.trial
            * j.current
            / Rloc
        )

        self.assertEqual(sp.expand(result.rhs - expected_rhs), 0)
        self.assertEqual(sp.expand(result.amat[psi] - expected_psi), 0)
        self.assertEqual(sp.expand(result.amat[j] - expected_j), 0)
        self.assertEqual(sp.expand(result.amat[T] - expected_T), 0)

    def test_constraint_uses_newton_residual_sign(self):
        psi = field("psi_constraint_test")
        j = field("j_constraint_test")
        v = test_function("v_constraint_test")
        Rloc = coefficient("R_constraint_test")
        C = (dR(v) * dR(psi.current) + v * j.current) / Rloc

        result = ConstraintEquation("current_definition_test", v, C).linearize(
            fields=(psi, j)
        )

        self.assertEqual(result.rhs, -C)
        self.assertEqual(result.kind, "static")
        self.assertEqual(result.amat[psi], dR(v) * dR(psi.trial) / Rloc)
        self.assertEqual(result.amat[j], v * j.trial / Rloc)



class FortranPrinterTest(unittest.TestCase):
    def test_model199_field_and_derivative_names(self):
        self.assertEqual(fortran(dR(j.current)), "zj0_x")
        self.assertEqual(fortran(dphi(u.trial)), "u_p")
        self.assertEqual(fortran(eta(psi.current)), "eta_T")

    def test_previous_increment_can_use_assembly_array_name(self):
        self.assertEqual(
            fortran(
                psi.previous_increment,
                previous_names={"psi": "delta_g(mp,1,ms,mt)"},
            ),
            "delta_g(mp,1,ms,mt)",
        )


if __name__ == "__main__":
    unittest.main()
