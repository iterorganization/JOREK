"""Dependency declarations for externally evaluated functions."""

from dataclasses import dataclass, field as dataclass_field
from enum import Enum
from typing import Dict, Mapping, Tuple, Union

import sympy as sp

from .exceptions import InvalidDeclarationError
from .symbols import validate_name


class LinearizationPolicy(str, Enum):
    ACTIVE = "active"
    FROZEN = "frozen"
    PIECEWISE_ACTIVE = "piecewise_active"


_FUNCTION_REGISTRY: Dict[object, "ExternalFunction"] = {}
_NAME_REGISTRY: Dict[str, "ExternalFunction"] = {}


@dataclass(frozen=True)
class ExternalFunction:
    """Declaration of a value supplied outside the symbolic equation layer."""

    name: str
    arguments: Tuple[str, ...]
    derivatives: Mapping[str, str]
    policy: LinearizationPolicy = LinearizationPolicy.ACTIVE
    fortran_name: str = ""
    _sympy_function: object = dataclass_field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        validate_name(self.name, "external-function")
        arguments = tuple(self.arguments)
        if len(set(arguments)) != len(arguments):
            raise InvalidDeclarationError(
                "External function {} has duplicate arguments".format(self.name)
            )
        for argument in arguments:
            validate_name(argument, "external-function argument")

        derivatives = dict(self.derivatives)
        unknown = set(derivatives) - set(arguments)
        if unknown:
            raise InvalidDeclarationError(
                "Derivatives for {} refer to unknown arguments: {}".format(
                    self.name, ", ".join(sorted(unknown))
                )
            )
        for derivative_name in derivatives.values():
            validate_name(derivative_name, "external derivative")

        policy = LinearizationPolicy(self.policy)
        fortran_name = self.fortran_name or self.name
        validate_name(fortran_name, "Fortran external-function")

        object.__setattr__(self, "arguments", arguments)
        object.__setattr__(self, "derivatives", derivatives)
        object.__setattr__(self, "policy", policy)
        object.__setattr__(self, "fortran_name", fortran_name)
        object.__setattr__(self, "_sympy_function", sp.Function(self.name))

    def __call__(self, *values):
        if len(values) != len(self.arguments):
            raise InvalidDeclarationError(
                "{} expects {} arguments ({}), received {}".format(
                    self.name,
                    len(self.arguments),
                    ", ".join(self.arguments),
                    len(values),
                )
            )
        return self._sympy_function(*(sp.sympify(value) for value in values))

    def derivative_call(self, argument: str, values):
        """Return the symbolic value of one supplied partial derivative."""

        return sp.Function(self.derivatives[argument])(*values)


def external_function(
    name: str,
    *,
    arguments: Tuple[str, ...],
    derivatives: Mapping[str, str],
    policy: Union[str, LinearizationPolicy] = LinearizationPolicy.ACTIVE,
    fortran_name: str = ""
) -> ExternalFunction:
    """Declare and register an externally evaluated symbolic function."""

    candidate = ExternalFunction(
        name=name,
        arguments=arguments,
        derivatives=derivatives,
        policy=LinearizationPolicy(policy),
        fortran_name=fortran_name,
    )
    existing = _NAME_REGISTRY.get(name)
    if existing is not None:
        if existing != candidate:
            raise InvalidDeclarationError(
                "External function {!r} was redeclared inconsistently".format(name)
            )
        return existing

    _NAME_REGISTRY[name] = candidate
    _FUNCTION_REGISTRY[candidate._sympy_function] = candidate
    return candidate


def definition_for_call(expression):
    """Return the registered declaration for a SymPy function call."""

    return _FUNCTION_REGISTRY.get(expression.func)

