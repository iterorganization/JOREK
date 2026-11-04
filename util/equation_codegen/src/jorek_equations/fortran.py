"""Fortran-style rendering of the symbolic DSL expressions."""

import sympy as sp
from sympy.printing.str import StrPrinter

from .external import definition_for_call
from .operators import R, Z, phi, s, t, PoloidalBracket, SpatialDerivative, expand_derivatives
from .symbols import FieldRole, FieldValue, Frozen, TestFunction


DEFAULT_FIELD_NAMES = {
    "psi": ("ps0", "psi"),
    "u": ("u0", "u"),
    "j": ("zj0", "zj"),
    "omega": ("w0", "w"),
    "rho": ("r0", "rho"),
    "T": ("T0", "T"),
}

COORDINATE_NAMES = {R: "BigR", Z: "Z", phi: "phi"}
DERIVATIVE_SUFFIXES = {R: "x", Z: "y", phi: "p", s: "s", t: "t"}


class FortranPrinter(StrPrinter):
    """Render the current model-199 symbolic naming convention."""

    def __init__(
        self, *, field_names=None, coordinate_names=None, previous_names=None, **settings
    ):
        super().__init__(settings)
        self.field_names = dict(DEFAULT_FIELD_NAMES)
        if field_names:
            self.field_names.update(field_names)
        self.coordinate_names = dict(COORDINATE_NAMES)
        if coordinate_names:
            self.coordinate_names.update(coordinate_names)
        self.previous_names = dict(previous_names or {})

    def _print(self, expr, **kwargs):
        # DSL atoms provide _sympystr methods for normal SymPy output. Intercept
        # them here so the Fortran backend can apply its own naming convention.
        if isinstance(expr, SpatialDerivative):
            return self._print_SpatialDerivative(expr)
        if isinstance(expr, PoloidalBracket):
            return self._print_PoloidalBracket(expr)
        if isinstance(expr, FieldValue):
            return self._print_FieldValue(expr)
        if isinstance(expr, Frozen):
            return self._print_Frozen(expr)
        if isinstance(expr, TestFunction):
            return self._print_TestFunction(expr)
        if isinstance(expr, sp.Function):
            return self._print_Function(expr)
        return super()._print(expr, **kwargs)

    def _print_Symbol(self, expr):
        return self.coordinate_names.get(expr, super()._print_Symbol(expr))

    def _print_TestFunction(self, expr):
        return expr.name

    def _print_Frozen(self, expr):
        # Frozen values are evaluated by the host before assembly.  Do not
        # expose the Python ``freeze(Abs(...))`` wrapper in generated Fortran.
        return self._print(expr.expression)

    def _print_Abs(self, expr):
        # Model-199 stores abs(rho) and abs(T) in the current-state arrays;
        # their DSL derivatives are intentionally disabled.
        return self._print(expr.args[0])

    def _print_FieldValue(self, expr):
        current, trial = self.field_names.get(
            expr.field_name, (expr.field_name + "0", expr.field_name)
        )
        if expr.role is FieldRole.CURRENT:
            return current
        if expr.role is FieldRole.TRIAL:
            return trial
        if expr.role is FieldRole.PREVIOUS_DELTA:
            if expr.field_name in self.previous_names:
                return self.previous_names[expr.field_name]
            return "delta_{}_prev".format(trial)
        if expr.role is FieldRole.DELTA:
            return "delta_{}".format(trial)
        return trial

    def _print_SpatialDerivative(self, expr):
        base, suffix = self._derivative_parts(expr)
        return self._print(base) + "_" + suffix

    def _print_PoloidalBracket(self, expr):
        return "bracket({}, {})".format(
            self._print(expr.left), self._print(expr.right)
        )

    def _derivative_parts(self, expr):
        if isinstance(expr, SpatialDerivative):
            base, suffix = self._derivative_parts(expr.expression)
            suffix += DERIVATIVE_SUFFIXES[expr.coordinate]
            # JOREK writes mixed derivatives with poloidal derivatives first
            # and toroidal derivatives last (``u_xpp``, not ``u_ppx``).
            order = {"x": 0, "y": 1, "s": 2, "t": 3, "p": 4}
            suffix = "".join(sorted(suffix, key=lambda item: order[item]))
            return base, suffix
        return expr, ""

    def _print_Function(self, expr):
        external = definition_for_call(expr)
        if external is not None:
            return external.fortran_name
        # Supplied derivative functions are represented as SymPy calls in the
        # symbolic tree, but JOREK evaluates them into scalar variables.
        for definition in self._external_definitions():
            if expr.func.__name__ in definition.derivatives.values():
                return expr.func.__name__
        return super()._print_Function(expr)

    def _external_definitions(self):
        # The registry is intentionally queried lazily to keep this printer's
        # public API independent of registry internals.
        from .external import _NAME_REGISTRY

        return tuple(_NAME_REGISTRY.values())


def fortran(expression, **settings) -> str:
    """Return a JOREK-style Fortran expression string."""

    return FortranPrinter(**settings).doprint(expand_derivatives(expression))
