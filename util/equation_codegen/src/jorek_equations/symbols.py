"""Core symbolic objects used by the JOREK equation DSL."""

from dataclasses import dataclass
from enum import Enum
import re
from typing import Optional, Union

import sympy as sp

from .exceptions import InvalidDeclarationError


_VALID_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def validate_name(name: str, kind: str) -> None:
    if not isinstance(name, str) or not _VALID_NAME.match(name):
        raise InvalidDeclarationError(
            "Invalid {} name {!r}; expected an identifier".format(kind, name)
        )


class FieldRole(str, Enum):
    """Role played by a field value during linearization."""

    CURRENT = "current"
    ABSTRACT = "abstract"
    DELTA = "delta"
    PREVIOUS_DELTA = "previous_delta"
    TRIAL = "trial"


class FieldValue(sp.Expr):
    """One physical field in a specific linearization role."""

    is_commutative = True
    is_Atom = True

    def __new__(cls, field_name: str, role: Union[str, FieldRole]):
        validate_name(field_name, "field")
        role_value = FieldRole(role).value
        return sp.Expr.__new__(
            cls, sp.Symbol(field_name), sp.Symbol(role_value)
        )

    @property
    def field_name(self) -> str:
        return str(self.args[0])

    @property
    def role(self) -> FieldRole:
        return FieldRole(str(self.args[1]))

    def _sympystr(self, printer) -> str:
        if self.role is FieldRole.ABSTRACT:
            return self.field_name
        if self.role is FieldRole.CURRENT:
            return "{}0".format(self.field_name)
        if self.role is FieldRole.DELTA:
            return "delta_{}".format(self.field_name)
        if self.role is FieldRole.PREVIOUS_DELTA:
            return "delta_{}_prev".format(self.field_name)
        return "{}_trial".format(self.field_name)

    def _latex(self, printer) -> str:
        name = printer._print(sp.Symbol(self.field_name))
        if self.role is FieldRole.ABSTRACT:
            return name
        if self.role is FieldRole.CURRENT:
            return "{}_0".format(name)
        if self.role is FieldRole.DELTA:
            return r"\delta {}".format(name)
        if self.role is FieldRole.PREVIOUS_DELTA:
            return r"\delta {}^{{n-1}}".format(name)
        return r"\varphi_{{{}}}".format(name)


@dataclass(frozen=True)
class Field:
    """Declaration of one evolved physical field."""

    name: str
    fortran_current: Optional[str] = None
    fortran_trial: Optional[str] = None

    def __post_init__(self) -> None:
        validate_name(self.name, "field")
        if self.fortran_current is not None:
            validate_name(self.fortran_current, "Fortran current-value")
        if self.fortran_trial is not None:
            validate_name(self.fortran_trial, "Fortran trial-value")

    @property
    def current(self) -> FieldValue:
        return FieldValue(self.name, FieldRole.CURRENT)

    @property
    def symbol(self) -> FieldValue:
        """Return the abstract symbol used in equation definitions."""

        return FieldValue(self.name, FieldRole.ABSTRACT)

    @property
    def increment(self) -> FieldValue:
        return FieldValue(self.name, FieldRole.DELTA)

    @property
    def previous_increment(self) -> FieldValue:
        return FieldValue(self.name, FieldRole.PREVIOUS_DELTA)

    @property
    def trial(self) -> FieldValue:
        return FieldValue(self.name, FieldRole.TRIAL)

    def value(self, role: Union[str, FieldRole]) -> FieldValue:
        return FieldValue(self.name, role)

    def _sympy_(self):
        """Allow a bare field declaration in SymPy expressions."""

        return self.symbol

    # Delegate arithmetic to the abstract symbolic field, so equations can be
    # written as ``v * eta(T) * j / R``.
    def __add__(self, other):
        return self.symbol + other

    def __radd__(self, other):
        return other + self.symbol

    def __sub__(self, other):
        return self.symbol - other

    def __rsub__(self, other):
        return other - self.symbol

    def __mul__(self, other):
        return self.symbol * other

    def __rmul__(self, other):
        return other * self.symbol

    def __truediv__(self, other):
        return self.symbol / other

    def __rtruediv__(self, other):
        return other / self.symbol

    def __pow__(self, other):
        return self.symbol ** other

    def __neg__(self):
        return -self.symbol


class TestFunction(sp.Expr):
    """Symbolic Galerkin test function."""

    is_commutative = True
    is_Atom = True

    def __new__(cls, name: str):
        validate_name(name, "test-function")
        return sp.Expr.__new__(cls, sp.Symbol(name))

    @property
    def name(self) -> str:
        return str(self.args[0])

    def _sympystr(self, printer) -> str:
        return self.name

    def _latex(self, printer) -> str:
        return printer._print(sp.Symbol(self.name))


class Frozen(sp.Expr):
    """Expression evaluated at the current state but excluded from variation."""

    is_commutative = True

    def __new__(cls, expression):
        return sp.Expr.__new__(cls, sp.sympify(expression))

    @property
    def expression(self):
        return self.args[0]

    def _sympystr(self, printer) -> str:
        return "freeze({})".format(printer._print(self.expression))


def field(
    name: str,
    *,
    fortran_current: Optional[str] = None,
    fortran_trial: Optional[str] = None
) -> Field:
    """Declare an evolved physical field."""

    return Field(name, fortran_current, fortran_trial)


def test_function(name: str = "v") -> TestFunction:
    """Declare a symbolic weak-form test function."""

    return TestFunction(name)


def coefficient(name: str, **assumptions) -> sp.Symbol:
    """Declare a state-independent scalar coefficient."""

    validate_name(name, "coefficient")
    return sp.Symbol(name, **assumptions)


def delta(value: Field) -> FieldValue:
    return value.increment


def previous_delta(value: Field) -> FieldValue:
    return value.previous_increment


def trial(value: Field) -> FieldValue:
    return value.trial


def freeze(expression) -> Frozen:
    """Mark a value as frozen during directional linearization."""

    return Frozen(expression)


def frozen_abs(expression) -> Frozen:
    """Return an absolute value evaluated at the current state only."""

    return freeze(sp.Abs(sp.sympify(expression)))
