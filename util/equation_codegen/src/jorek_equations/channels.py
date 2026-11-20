"""Classification by toroidal derivatives on test and trial functions."""

from dataclasses import dataclass

import sympy as sp

from .operators import SpatialDerivative, expand_derivatives, phi
from .symbols import FieldRole, FieldValue, TestFunction


@dataclass(frozen=True)
class ToroidalChannels:
    """The four JOREK FFT assembly channels for one expression."""

    p: sp.Expr = sp.S.Zero
    n: sp.Expr = sp.S.Zero
    k: sp.Expr = sp.S.Zero
    kn: sp.Expr = sp.S.Zero


def _term_channel(term):
    test_order = 0
    trial_order = 0
    for item in term.atoms(SpatialDerivative):
        if item.coordinate != phi:
            continue
        if item.expression.has(TestFunction):
            test_order += 1
        if any(
            value.role is FieldRole.TRIAL
            for value in item.expression.atoms(FieldValue)
        ):
            trial_order += 1

    if test_order > 1 or trial_order > 1:
        raise ValueError(
            "Unsupported toroidal derivative order in term {}".format(term)
        )
    if test_order and trial_order:
        return "kn"
    if test_order:
        return "k"
    if trial_order:
        return "n"
    return "p"


def split_toroidal_channels(expression) -> ToroidalChannels:
    """Split an expanded sum into p, n, k, and kn contributions."""

    # Apply the product rule before classifying a term.  For example,
    # ``dphi(T * delta_rho)`` must become ``T_phi*delta_rho`` (p channel)
    # plus ``T*dphi(delta_rho)`` (n channel).  Classifying the unevaluated
    # derivative would incorrectly send both contributions to ``n``.
    expression = expand_derivatives(expression)
    groups = {"p": [], "n": [], "k": [], "kn": []}
    for term in sp.Add.make_args(sp.expand(expression)):
        groups[_term_channel(term)].append(term)
    return ToroidalChannels(
        **{name: sp.Add(*terms) for name, terms in groups.items()}
    )
