"""Exceptions raised by the JOREK equation DSL."""


class EquationDslError(Exception):
    """Base class for errors produced by the equation DSL."""


class InvalidDeclarationError(EquationDslError, ValueError):
    """A symbolic object was declared inconsistently."""


class MissingDerivativeError(EquationDslError):
    """An active external dependency has no supplied derivative."""


class UnsupportedExpressionError(EquationDslError):
    """An expression cannot be linearized safely."""

