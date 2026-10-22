"""Directional linearization of JOREK DSL expressions."""

from functools import lru_cache

import sympy as sp
from sympy.core.function import AppliedUndef, ArgumentIndexError

from .exceptions import MissingDerivativeError, UnsupportedExpressionError
from .external import LinearizationPolicy, definition_for_call
from .operators import PoloidalBracket, SpatialDerivative, bracket, derivative
from .symbols import Field, FieldRole, FieldValue, Frozen, TestFunction


def _is_zero(expression) -> bool:
    return expression == 0


def resolve_current(expression):
    """Resolve abstract equation fields to current-state values."""

    expression = sp.sympify(expression)
    replacements = {
        value: FieldValue(value.field_name, FieldRole.CURRENT)
        for value in expression.atoms(FieldValue)
        if value.role is FieldRole.ABSTRACT
    }
    return expression.xreplace(replacements)


def variation(expression, wrt: Field, *, direction=FieldRole.DELTA):
    """Compute the directional variation of an expression for one field."""

    if not isinstance(wrt, Field):
        raise TypeError("wrt must be a Field declaration")
    direction = FieldRole(direction)
    if direction is FieldRole.CURRENT:
        raise ValueError("A variation direction cannot have role 'current'")
    # Equation definitions contain abstract fields. At the tangent stage,
    # non-varied factors are evaluated at the current state.
    expression = resolve_current(expression)

    @lru_cache(maxsize=None)
    def visit(node):
        if isinstance(node, Frozen):
            return sp.S.Zero

        if isinstance(node, FieldValue):
            if node.role is FieldRole.CURRENT and node.field_name == wrt.name:
                return wrt.value(direction)
            return sp.S.Zero

        if isinstance(node, (TestFunction, sp.Symbol, sp.Number)):
            return sp.S.Zero

        if isinstance(node, SpatialDerivative):
            inner = visit(node.expression)
            if _is_zero(inner):
                return sp.S.Zero
            return derivative(inner, node.coordinate)

        if isinstance(node, PoloidalBracket):
            left_variation = visit(node.left)
            right_variation = visit(node.right)
            result = sp.S.Zero
            if not _is_zero(left_variation):
                result += bracket(left_variation, node.right)
            if not _is_zero(right_variation):
                result += bracket(node.left, right_variation)
            return result

        external = definition_for_call(node)
        if external is not None:
            if external.policy is LinearizationPolicy.FROZEN:
                return sp.S.Zero
            result = sp.S.Zero
            for argument_name, argument_value in zip(external.arguments, node.args):
                argument_variation = visit(argument_value)
                if _is_zero(argument_variation):
                    continue
                if argument_name not in external.derivatives:
                    raise MissingDerivativeError(
                        "External function {} is active in argument {!r}, but no "
                        "derivative was supplied".format(external.name, argument_name)
                    )
                result += (
                    external.derivative_call(argument_name, node.args)
                    * argument_variation
                )
            return result

        if isinstance(node, AppliedUndef):
            raise UnsupportedExpressionError(
                "Undefined function {} has no external-function declaration".format(
                    node.func
                )
            )

        if isinstance(node, sp.Add):
            return sp.Add(*(visit(argument) for argument in node.args))

        if isinstance(node, sp.Mul):
            terms = []
            for index, argument in enumerate(node.args):
                argument_variation = visit(argument)
                if _is_zero(argument_variation):
                    continue
                factors = list(node.args)
                factors[index] = argument_variation
                terms.append(sp.Mul(*factors))
            return sp.Add(*terms)

        if isinstance(node, sp.Pow):
            base, exponent = node.args
            base_variation = visit(base)
            exponent_variation = visit(exponent)
            result = sp.S.Zero
            if not _is_zero(base_variation):
                result += exponent * base ** (exponent - 1) * base_variation
            if not _is_zero(exponent_variation):
                result += node * exponent_variation * sp.log(base)
            return result

        if isinstance(node, sp.Function):
            result = sp.S.Zero
            for index, argument in enumerate(node.args, start=1):
                argument_variation = visit(argument)
                if _is_zero(argument_variation):
                    continue
                try:
                    partial = node.fdiff(index)
                except (ArgumentIndexError, NotImplementedError) as exc:
                    raise UnsupportedExpressionError(
                        "Cannot differentiate function {} safely".format(node.func)
                    ) from exc
                result += partial * argument_variation
            return result

        if not node.args:
            return sp.S.Zero

        raise UnsupportedExpressionError(
            "Unsupported expression node {} in {}".format(type(node).__name__, node)
        )

    return sp.expand(visit(expression))


def jacobian_blocks(expression, fields):
    """Return one trial-function directional derivative per field."""

    return {
        value: variation(expression, value, direction=FieldRole.TRIAL)
        for value in fields
    }
