"""Weak equation containers and JOREK sign transformations."""

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping

import sympy as sp

from .linearization import resolve_current, variation
from .symbols import Field, FieldRole, TestFunction


@dataclass(frozen=True)
class LinearizedEquation:
    """Symbolic local RHS and AMAT blocks for one weak equation."""

    name: str
    rhs: sp.Expr
    amat: Dict[Field, sp.Expr]
    kind: str = "evolution"


@dataclass(frozen=True)
class EvolutionEquation:
    """Weak evolution equation dA/dt = B."""

    name: str
    test: TestFunction
    A: sp.Expr
    B: sp.Expr
    kind: str = "evolution"
    # Optional field-specific residual expressions used when a legacy model
    # defines a Jacobian tangent that intentionally differs from dB/dq.
    amat_variation_overrides: Mapping[Field, sp.Expr] = None

    def __post_init__(self) -> None:
        if not isinstance(self.test, TestFunction):
            raise TypeError("test must be a TestFunction")
        object.__setattr__(self, "A", sp.sympify(self.A))
        object.__setattr__(self, "B", sp.sympify(self.B))
        overrides = self.amat_variation_overrides or {}
        object.__setattr__(
            self,
            "amat_variation_overrides",
            {field: sp.sympify(expression) for field, expression in overrides.items()},
        )

    def linearize(self, *, fields: Iterable[Field], timestep, theta, zeta):
        """Apply the model-199 evolution-equation sign convention."""

        fields = tuple(fields)
        timestep = sp.sympify(timestep)
        theta = sp.sympify(theta)
        zeta = sp.sympify(zeta)

        history = sp.Add(
            *(
                variation(self.A, value, direction=FieldRole.PREVIOUS_DELTA)
                for value in fields
            )
        )
        rhs = sp.expand(timestep * resolve_current(self.B) + zeta * history)
        amat = {
            value: sp.expand(
                (1 + zeta)
                * variation(self.A, value, direction=FieldRole.TRIAL)
                - theta
                * timestep
                * variation(
                    self.amat_variation_overrides.get(value, self.B),
                    value,
                    direction=FieldRole.TRIAL,
                )
            )
            for value in fields
        }
        return LinearizedEquation(self.name, rhs, amat, self.kind)


@dataclass(frozen=True)
class ConstraintEquation:
    """Weak algebraic constraint C(q) = 0."""

    name: str
    test: TestFunction
    C: sp.Expr
    kind: str = "static"

    def __post_init__(self) -> None:
        if not isinstance(self.test, TestFunction):
            raise TypeError("test must be a TestFunction")
        object.__setattr__(self, "C", sp.sympify(self.C))

    def linearize(self, *, fields: Iterable[Field]):
        """Return RHS=-C and AMAT=D(C) for the Newton correction."""

        fields = tuple(fields)
        rhs = -resolve_current(self.C)
        amat = {
            value: variation(self.C, value, direction=FieldRole.TRIAL)
            for value in fields
        }
        return LinearizedEquation(self.name, rhs, amat, self.kind)
