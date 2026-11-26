"""Model 199: read the element routine, generate the weak form, compare."""

import unittest
from pathlib import Path

import sympy as sp

from jorek_equations import (
    coefficient,
    dR,
    dZ,
    dphi,
    expand_brackets,
    split_toroidal_channels,
    variation,
)
from jorek_equations.model199 import (
    FIELDS,
    T,
    F0,
    eta,
    induction_equation_1,
    j,
    omega,
    psi,
    rho,
    u,
)
from jorek_equations.fortran_source import (
    EQUATION_1_ASSIGNMENTS,
    FortranComparisonError,
    compare_assignment_maps,
    compare_model199_equation1,
    extract_fortran_assignments,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
FORTRAN_FILE = REPOSITORY_ROOT / "models/model199/mod_elt_matrix_fft.f90"


class SourceComparisonTest(unittest.TestCase):
    def test_extracts_all_equation_1_assignments(self):
        assignments = extract_fortran_assignments(
            FORTRAN_FILE, EQUATION_1_ASSIGNMENTS
        )
        self.assertEqual(set(assignments), set(EQUATION_1_ASSIGNMENTS))
        self.assertIn("eta_T", assignments["rhs_ij_1"])
        self.assertIn("deta_dT", assignments["amat_16"])

    def test_real_model199_equation_1_matches(self):
        compare_model199_equation1(FORTRAN_FILE)

    def test_toroidal_trial_term_is_in_n_channel(self):
        result = induction_equation_1().linearize(
            fields=FIELDS,
            timestep=coefficient("tstep"),
            theta=coefficient("theta"),
            zeta=coefficient("zeta"),
        )
        channels = split_toroidal_channels(result.amat[u])
        self.assertNotEqual(channels.p, 0)
        self.assertNotEqual(channels.n, 0)
        self.assertEqual(channels.k, 0)
        self.assertEqual(channels.kn, 0)

    def test_mismatch_raises_and_reports_difference(self):
        with self.assertRaises(FortranComparisonError) as raised:
            compare_assignment_maps(
                {"amat_test": "x + y"},
                {"amat_test": "x - y"},
            )
        message = str(raised.exception)
        self.assertIn("amat_test", message)
        self.assertIn("source-only:", message)
        self.assertIn("generated-only:", message)
        self.assertIn("      + y", message)
        self.assertIn("      - y", message)


if __name__ == "__main__":
    unittest.main()

class Model199Equation1Test(unittest.TestCase):
    def setUp(self):
        self.timestep = coefficient("tstep_model199_eq1_test")
        self.theta = coefficient("theta_model199_eq1_test")
        self.zeta = coefficient("zeta_model199_eq1_test")
        self.xjac = coefficient("xjac")
        self.R = sp.Symbol("R")

    def test_rhs_matches_fortran_equation_1_terms(self):
        equation = induction_equation_1()
        result = equation.linearize(
            fields=(psi, u, j, T),
            timestep=self.timestep,
            theta=self.theta,
            zeta=self.zeta,
        )
        v = equation.test
        current_source = sp.Symbol("current_source")
        eta_num = sp.Symbol("eta_num")
        eps_cyl = sp.Symbol("eps_cyl")

        expected_B = (
            v * eta(T.current) * (j.current - current_source) / self.R * self.xjac
            + v * self.xjac * (
                dR(psi.current) * dZ(u.current)
                - dZ(psi.current) * dR(u.current)
            )
            - v * eps_cyl * F0 / self.R * dphi(u.current) * self.xjac
            + eta_num * self.xjac * (
                dR(v) * dR(j.current) + dZ(v) * dZ(j.current)
            )
        )
        expected = self.timestep * expected_B + self.zeta * (
            v * psi.previous_increment / self.R * self.xjac
        )
        self.assertEqual(sp.expand(expand_brackets(result.rhs) - expected), 0)

    def test_amat_blocks_match_fortran_equation_1_terms(self):
        equation = induction_equation_1()
        result = equation.linearize(
            fields=(psi, u, j, T),
            timestep=self.timestep,
            theta=self.theta,
            zeta=self.zeta,
        )
        v = equation.test
        eta_num = sp.Symbol("eta_num")
        eps_cyl = sp.Symbol("eps_cyl")
        current_source = sp.Symbol("current_source")

        expected_psi = (
            (1 + self.zeta) * v * psi.trial / self.R * self.xjac
            - self.theta * self.timestep * v * self.xjac * (
                dR(psi.trial) * dZ(u.current)
                - dZ(psi.trial) * dR(u.current)
            )
        )
        expected_u = (
            -self.theta * self.timestep * v * self.xjac * (
                dR(psi.current) * dZ(u.trial)
                - dZ(psi.current) * dR(u.trial)
            )
            + self.theta * self.timestep * v * eps_cyl * F0 / self.R
            * dphi(u.trial) * self.xjac
        )
        expected_j = -self.theta * self.timestep * (
            eta_num * self.xjac * (
                dR(v) * dR(j.trial) + dZ(v) * dZ(j.trial)
            )
            + eta(T.current) * v * j.trial / self.R * self.xjac
        )
        expected_T = -self.theta * self.timestep * (
            sp.Function("deta_dT")(T.current) * T.trial * v
            * (j.current - current_source) / self.R * self.xjac
        )

        self.assertEqual(
            sp.expand(expand_brackets(result.amat[psi]) - expected_psi), 0
        )
        self.assertEqual(
            sp.expand(expand_brackets(result.amat[u]) - expected_u), 0
        )
        self.assertEqual(sp.expand(result.amat[j] - expected_j), 0)
        self.assertEqual(sp.expand(result.amat[T] - expected_T), 0)

    def test_unrelated_model_fields_have_zero_blocks(self):
        equation = induction_equation_1()
        result = equation.linearize(
            fields=FIELDS,
            timestep=self.timestep,
            theta=self.theta,
            zeta=self.zeta,
        )
        self.assertEqual(result.amat[omega], 0)
        self.assertEqual(result.amat[rho], 0)


if __name__ == "__main__":
    unittest.main()


if __name__ == "__main__":
    unittest.main()
