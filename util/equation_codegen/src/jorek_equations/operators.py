"""Structured spatial operators for the equation DSL."""

from typing import Sequence, Tuple

import sympy as sp

from .external import definition_for_call

from .exceptions import InvalidDeclarationError
from .symbols import Field, FieldValue, Frozen, TestFunction


R = sp.Symbol("R")
Z = sp.Symbol("Z")
phi = sp.Symbol("phi")
s = sp.Symbol("s")
t = sp.Symbol("t")


class SpatialDerivative(sp.Expr):
    """An unevaluated derivative with respect to a physical coordinate."""

    is_commutative = True

    def __new__(cls, expression, coordinate):
        expression = sp.sympify(expression)
        coordinate = sp.sympify(coordinate)
        if coordinate not in (R, Z, phi, s, t):
            raise InvalidDeclarationError(
                "Unsupported physical coordinate {!r}".format(coordinate)
            )
        return sp.Expr.__new__(cls, expression, coordinate)

    @property
    def expression(self):
        return self.args[0]

    @property
    def coordinate(self):
        return self.args[1]

    def _sympystr(self, printer) -> str:
        return "d{}({})".format(
            printer._print(self.coordinate), printer._print(self.expression)
        )

    def _latex(self, printer) -> str:
        return r"\frac{{\partial {}}}{{\partial {}}}".format(
            printer._print(self.expression), printer._print(self.coordinate)
        )


def derivative(expression, coordinate) -> SpatialDerivative:
    return SpatialDerivative(expression, coordinate)


def dR(expression) -> SpatialDerivative:
    return derivative(expression, R)


def dZ(expression) -> SpatialDerivative:
    return derivative(expression, Z)


def dphi(expression) -> SpatialDerivative:
    return derivative(expression, phi)


def ds(expression) -> SpatialDerivative:
    return derivative(expression, s)


def dt(expression) -> SpatialDerivative:
    return derivative(expression, t)


def grad(expression) -> Tuple[SpatialDerivative, SpatialDerivative]:
    """Return the two poloidal gradient components."""

    return dR(expression), dZ(expression)


def dot(left: Sequence, right: Sequence):
    if len(left) != len(right):
        raise InvalidDeclarationError("Dot-product operands have different sizes")
    return sp.Add(*(sp.sympify(a) * sp.sympify(b) for a, b in zip(left, right)))


class PoloidalBracket(sp.Expr):
    """Structured physical-coordinate bracket [a,b]."""

    is_commutative = True

    def __new__(cls, left, right):
        return sp.Expr.__new__(cls, sp.sympify(left), sp.sympify(right))

    @property
    def left(self):
        return self.args[0]

    @property
    def right(self):
        return self.args[1]

    def _sympystr(self, printer) -> str:
        return "bracket({}, {})".format(
            printer._print(self.left), printer._print(self.right)
        )


def bracket(left, right):
    """Return the structured poloidal bracket [a,b]."""

    return PoloidalBracket(left, right)


def poiss_bracket_st(left, right):
    """Return ``left_s right_t - left_t right_s`` in element coordinates."""

    return ds(left) * dt(right) - dt(left) * ds(right)


def poiss_bracket(left, right):
    """Return ``left_R right_Z - left_Z right_R`` in physical coordinates."""

    return dR(left) * dZ(right) - dZ(left) * dR(right)


def laplacian(value):
    """Return the axisymmetric (R,Z) Laplacian of ``value``."""

    return dR(dR(value)) + dZ(dZ(value)) + dR(value) / R


def expand_brackets(expression):
    """Expand physical brackets into R/Z derivatives."""

    expression = sp.sympify(expression)
    replacements = {
        item: dR(item.left) * dZ(item.right) - dZ(item.left) * dR(item.right)
        for item in expression.atoms(PoloidalBracket)
    }
    return expression.xreplace(replacements)


def expand_derivatives(expression):
    """Apply the product rule to structured spatial derivatives."""

    expression = sp.sympify(expression)

    def expand(node):
        if isinstance(node, (FieldValue, Frozen, TestFunction)):
            return node
        if isinstance(node, SpatialDerivative):
            inner = expand(node.expression)
            coordinate = node.coordinate
            if inner.is_Number:
                # A numeric factor (such as the 1/2 in the one-temperature
                # ion temperature Ti0 = T0/2) is constant in space.  Without
                # this the product rule below emits a spurious derivative of
                # that number.
                return sp.S.Zero
            if isinstance(inner, sp.Add):
                return sp.Add(*(expand(SpatialDerivative(term, coordinate))
                                for term in inner.args))
            if isinstance(inner, sp.Mul):
                terms = []
                for index, factor in enumerate(inner.args):
                    factors = list(inner.args)
                    factors[index] = expand(SpatialDerivative(factor, coordinate))
                    terms.append(sp.Mul(*factors))
                return sp.Add(*terms)
            if isinstance(inner, sp.Pow) and inner.exp.is_integer:
                return inner.exp * inner.base ** (inner.exp - 1) * expand(
                    SpatialDerivative(inner.base, coordinate)
                )
            external = definition_for_call(inner)
            if external is not None:
                return sp.Add(*(
                    external.derivative_call(argument_name, inner.args)
                    * expand(SpatialDerivative(argument_value, coordinate))
                    for argument_name, argument_value in zip(
                        external.arguments, inner.args
                    )
                    if argument_name in external.derivatives
                ))
            if inner in (R, Z, phi, s, t):
                return sp.S.One if inner == coordinate else sp.S.Zero
            return SpatialDerivative(inner, coordinate)
        if not node.args:
            return node
        return node.func(*[expand(argument) for argument in node.args])

    return sp.expand(expand(expression))


def lower_brackets_to_element(expression, jacobian):
    """Rewrite [a,b] as (a_s b_t-a_t b_s)/jacobian."""

    expression = sp.sympify(expression)
    jacobian = sp.sympify(jacobian)
    replacements = {
        item: (
            ds(item.left) * dt(item.right)
            - dt(item.left) * ds(item.right)
        )
        / jacobian
        for item in expression.atoms(PoloidalBracket)
    }
    return sp.expand(expression.xreplace(replacements))
