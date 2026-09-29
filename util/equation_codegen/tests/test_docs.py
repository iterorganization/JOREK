"""The model-600 documentation generated from model600.py."""

import ast
import sys
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "examples"))

from jorek_equations.latex_render import LatexPrinter  # noqa: E402
from jorek_equations.model600_docs import read_source, render_page  # noqa: E402
import model600_docs as tool  # noqa: E402


def latex(source, **options):
    return LatexPrinter({}, **options).latex(ast.parse(source, mode="eval").body)


class PrinterTest(unittest.TestCase):
    def test_quotient_keeps_its_factors_at_full_size(self):
        self.assertEqual(latex("2 * a / b * c / d"), r"\frac{2}{b\,d}\,a\,c")

    def test_element_bracket_absorbs_the_jacobian(self):
        self.assertEqual(latex("R * poiss_bracket_st(a, u) / xjac"), r"R\,[a,u]^{st}")
        self.assertEqual(latex("v * poiss_bracket_st(a, u) / dV"), r"\frac{1}{R}\,v\,[a,u]^{st}")

    def test_gradient_component_through_a_local(self):
        values = {"g": ast.parse("grad(psi)", mode="eval").body}
        self.assertEqual(latex("g[1]", values=values), r"\partial_{Z} \psi")

    def test_partial_alias_is_rendered_as_its_helper(self):
        self.assertEqual(latex("f(a)", aliases={"f": "dR"}), r"\partial_{R} a")

    def test_zip_comprehension_is_a_vector(self):
        self.assertEqual(
            latex("tuple(a - k * b for a, b in zip(grad(x), grad(y)))"),
            r"\nabla_{\mathrm{pol}} x - k\,\nabla_{\mathrm{pol}} y")


class Model600PageTest(unittest.TestCase):
    """weak_form.md is generated from model600.py and must follow it."""

    def test_every_equation_and_helper_is_documented(self):
        helpers, _, equations = read_source(tool.SOURCE.read_text(encoding="utf-8"))
        self.assertIn("neutral_density_equation_rhon", {e.name for e in equations})
        self.assertIn("_B_dot_grad", helpers)

    def test_page_is_up_to_date(self):
        text, new = tool.rendered()
        self.assertEqual(new, text, "run: python3 examples/model600_docs.py render")


if __name__ == "__main__":
    unittest.main()
