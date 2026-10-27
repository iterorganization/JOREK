"""Meld-friendly Markdown exports of model-600 source and generated terms."""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re

import sympy as sp

from .equations import LinearizedEquation
from .linearization import resolve_current, variation
from .fortran import fortran
from .model600 import (
    FIELDS,
    density_equation_rho,
    impurity_density_equation_rhoimp,
    neutral_density_equation_rhon,
    electron_energy_equation_Te,
    ion_energy_equation_Ti,
    total_energy_equation_T,
    parallel_velocity_equation_vpar,
    T,
    Te,
    current_constraint_equation_zj,
    induction_equation_1,
    j,
    momentum_equation_2,
    omega,
    psi,
    rho,
    u,
    vorticity_constraint_equation_w,
)
from .fortran_source import normalize_model600_text
from .operators import SpatialDerivative, expand_derivatives, phi
from .fortran_source import parse_fortran_expression
from .symbols import FieldRole, FieldValue, TestFunction, coefficient


ROWS = ("psi", "u", "zj", "w", "rho", "vpar", "rhoimp", "rhon", "ti", "te", "t")
FIELD_NAMES = {
    field: name
    for field, name in zip(FIELDS, ("psi", "u", "zj", "w", "rho", "T", "vpar", "Ti", "Te", "rhon", "rhoimp"))
}
ASSIGNMENT_RE = re.compile(
    r"(?P<lhs>(?:rhs_ij(?:_k)?|amat(?:_n|_k|_kn|_nn)?)\s*\(\s*"
    r"var_(?P<row>psi|u|zj|w|rho|rhoimp|rhon|vpar|Ti|Te|T)\b[^=]*?\))\s*=\s*(?P<rhs>.*)$",
    re.I,
)


@dataclass(frozen=True)
class SourceAssignment:
    lhs: str
    row: str
    line: int
    expression: str
    temperature_model: str = "both"


def _clean_lhs(lhs):
    return re.sub(r"\s+", "", lhs).lower()


def _read_source_assignments(path):
    """Read every relevant assignment, retaining repetitions and file order."""

    lines = Path(path).read_text(encoding="utf-8").splitlines()
    result = []
    temperature_model = "both"
    # Stack only conditionals controlled by ``with_TiTe``.  It is necessary
    # because model600 contains further ordinary IF blocks inside each
    # temperature branch, and a bare ``endif`` must not accidentally change
    # the active temperature convention.
    temperature_stack = []
    velocity_profile_branch = None
    index = 0
    while index < len(lines):
        raw = lines[index]
        code = raw.split("!", 1)[0]
        lowered = raw.lower()
        if re.search(r"\bif\s*\(\s*with_tite\s*\)\s*then", lowered):
            temperature_stack.append([temperature_model, 0])
            temperature_model = "two"
        elif temperature_stack and re.search(r"\bif\s*\([^)]*\)\s*then", lowered):
            temperature_stack[-1][1] += 1
        elif (
            temperature_stack
            and temperature_stack[-1][1] == 0
            and re.match(r"^\s*else\b(?!\s*if)", lowered)
        ):
            temperature_model = "single"
        elif temperature_stack and re.search(r"\bend\s*if\b", lowered):
            if temperature_stack[-1][1]:
                temperature_stack[-1][1] -= 1
            else:
                temperature_model = temperature_stack.pop()[0]
        # ``normalized_velocity_profile`` defaults to .true.; export that
        # branch and skip the ``else`` alternative.  These blocks contain no
        # nested conditionals, so a single flag is enough.
        if re.search(r"\bif\s*\(\s*normalized_velocity_profile\s*\)\s*then", lowered):
            velocity_profile_branch = True
            index += 1
            continue
        if velocity_profile_branch is not None:
            if re.match(r"^\s*else\b(?!\s*if)", lowered):
                velocity_profile_branch = False
                index += 1
                continue
            if re.search(r"\bend\s*if\b", lowered):
                velocity_profile_branch = None
                index += 1
                continue
            if not velocity_profile_branch:
                index += 1
                continue
        match = ASSIGNMENT_RE.search(code)
        if match is None:
            index += 1
            continue
        lhs = _clean_lhs(match.group("lhs"))
        pieces = []
        current = match.group("rhs")
        while True:
            current = current.strip()
            continued = current.endswith("&")
            pieces.append(current[:-1].strip() if continued else current)
            if not continued:
                break
            index += 1
            if index >= len(lines):
                raise ValueError("Unterminated assignment at line {}".format(index + 1))
            # JOREK often places a blank/comment line between continued
            # terms.  It is not the end of the Fortran statement.
            current = lines[index].split("!", 1)[0].strip()
            while not current:
                index += 1
                if index >= len(lines):
                    raise ValueError("Unterminated assignment at line {}".format(index + 1))
                current = lines[index].split("!", 1)[0].strip()
            if current.startswith("&"):
                current = current[1:]
        expression = " ".join(piece for piece in pieces if piece)
        # Repeated extension assignments have the form ``a = a + extension``.
        # Extension assignments commonly use either ``lhs = lhs + term`` or
        # ``lhs = lhs - term``.  Remove only the repeated lhs and preserve the
        # following sign so the extension remains a source term.
        expression = re.sub(
            r"^{}\s*".format(re.escape(lhs)), "",
            _clean_assignment_refs(expression), flags=re.I,
        )
        result.append(
            SourceAssignment(
                lhs, match.group("row").lower(), index + 2 - len(pieces),
                expression, temperature_model,
            )
        )
        index += 1
    return result


def _clean_assignment_refs(expression):
    """Canonicalize only assignment references, leaving term spelling intact."""

    return re.sub(
        r"(?:rhs_ij(?:_k)?|amat(?:_n|_k|_kn|_nn)?)\s*\([^)]*\)",
        lambda match: _clean_lhs(match.group(0)),
        expression,
        flags=re.I,
    )


def _split_top_level(expression):
    """Split an expression at ordered, outermost additive signs."""

    terms = []
    start = 0
    depth = 0
    for index, char in enumerate(expression):
        if char == "(":
            depth += 1
        elif char == ")":
            depth = max(0, depth - 1)
        elif char in "+-" and depth == 0 and index > start:
            previous = expression[index - 1]
            if previous not in "eEdD":
                terms.append(expression[start:index].strip())
                start = index
    tail = expression[start:].strip()
    if tail:
        terms.append(tail)
    # ``lhs = lhs + &`` continued onto a line that itself starts with ``+``
    # leaves a piece that is nothing but a sign.  It is not a term.
    return [term for term in terms if term.strip(" +-")]


def _normalized_text(text, *, single_temperature=False, two_temperature=False):
    text = re.sub(
        r"factor\(\s*var_[a-z0-9_]+\s*,\s*[0-9]+\s*\)", "1", text,
        flags=re.I,
    )
    text = re.sub(r"\bBigR_x\b", "1", text, flags=re.I)
    result = normalize_model600_text(
        text,
        single_temperature=single_temperature,
        two_temperature=two_temperature,
    )
    return re.sub(r"\bBigR_x\b", "1", result, flags=re.I)


def _parsed_terms(text, *, temperature_model=None):
    """Return expanded monomials, or the original text when unsupported."""

    variants = []
    variants_to_try = {
        "single": ((True, False),),
        "two": ((False, True),),
        None: ((False, False), (True, False), (False, True)),
    }[temperature_model]
    for single, two in variants_to_try:
        try:
            parsed = parse_fortran_expression(
                _normalized_text(text, single_temperature=single, two_temperature=two)
            )
            # ``normalize_fortran_text`` expands several element-routine
            # aliases (notably ``r0_x_hat``) after the textual normalization
            # above.  Those expansions can re-introduce ``BigR_x``.  The
            # report follows the comparison convention BigR_x=1, so apply it
            # to the parsed tree as well as to the input text.
            parsed = parsed.subs(sp.Symbol("BigR_x"), sp.Integer(1))
        except Exception:
            continue
        terms = tuple(sp.Add.make_args(sp.expand(parsed)))
        if terms not in variants:
            variants.append(terms)
    return variants


def _channel(term):
    test_order = 0
    trial_order = 0
    for derivative in term.atoms(SpatialDerivative):
        if derivative.coordinate != phi:
            continue
        if derivative.expression.has(TestFunction):
            test_order += 1
        if any(
            value.role is FieldRole.TRIAL
            for value in derivative.expression.atoms(FieldValue)
        ):
            trial_order += 1
    if test_order == 0 and trial_order == 0:
        return ""
    if test_order == 0 and trial_order == 1:
        return "_n"
    if test_order == 1 and trial_order == 0:
        return "_k"
    if test_order == 1 and trial_order == 1:
        return "_kn"
    if test_order == 0 and trial_order == 2:
        return "_nn"
    raise ValueError("Unsupported toroidal channel in {}".format(term))


def _add_linearized(
    pools, row, linearized: LinearizedEquation, *, previous_names=None,
    include_rhs=True, include_amat=True, include_fields=None,
):
    if include_rhs:
        for term in sp.Add.make_args(sp.expand(expand_derivatives(linearized.rhs))):
            if term == 0:
                continue
            channel = _channel(term)
            lhs = "rhs_ij{}(var_{})".format(channel, row)
            pools[_clean_lhs(lhs)].append(fortran(term, previous_names=previous_names))
    if include_amat:
        expressions = linearized.amat.items()
    else:
        expressions = ()
    if include_fields is not None:
        include_fields = set(include_fields)
        expressions = (
            (field, expression)
            for field, expression in expressions
            if field in include_fields
        )
    for field, expression in expressions:
        column = FIELD_NAMES[field]
        for term in sp.Add.make_args(sp.expand(expand_derivatives(expression))):
            if term == 0:
                continue
            channel = _channel(term)
            lhs = "amat{}(var_{},var_{})".format(channel, row, column)
            pools[_clean_lhs(lhs)].append(fortran(term, previous_names=previous_names))


def _raw_evolution_rhs(equation, fields, timestep, zeta):
    """Build RHS without SymPy's global additive expansion."""

    history = sp.Add(
        *(variation(equation.A, value, direction=FieldRole.PREVIOUS_DELTA) for value in fields)
    )
    terms = [
        timestep * resolve_current(term)
        for term in sp.Add.make_args(equation.B)
    ]
    terms.extend(zeta * term for term in sp.Add.make_args(history))
    return terms


def _add_raw_rhs(pools, row, expression, *, previous_names=None):
    for term in expression:
        if term == 0:
            continue
        channel = _channel(term)
        lhs = "rhs_ij{}(var_{})".format(channel, row)
        pools[_clean_lhs(lhs)].append(fortran(term, previous_names=previous_names))


def _generated_pools(
    rows=ROWS, *, include_neo=True, with_tite=None,
    st_form=True,
):
    """Generate implemented model-600 terms for the requested rows."""

    rows = set(rows)
    pools = defaultdict(list)
    timestep = coefficient("tstep")
    theta = coefficient("theta")
    zeta = coefficient("zeta")
    if "psi" in rows:
        # Both the residual and every tangent column follow the temperature
        # model of the section being exported.  The two branches differ in
        # more than naming: the one-temperature electron pressure is built
        # from Te0 = T0/2.
        psi_with_tite = bool(with_tite)
        equation = induction_equation_1(with_TiTe=psi_with_tite)
        result = equation.linearize(
            fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
        )
        _add_linearized(
            pools, "psi", result,
            previous_names={"psi": "delta_g(mp,var_psi,ms,mt)"},
        )
    if "rho" in rows:
        equation = density_equation_rho(with_TiTe=bool(with_tite))
        _add_linearized(
            pools, "rho",
            equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            ),
            previous_names={"rho": "delta_rho_g"},
        )
    if "vpar" in rows:
        equation = parallel_velocity_equation_vpar(
            with_TiTe=bool(with_tite),
            st_form=st_form,
        )
        _add_linearized(
            pools, "vpar",
            equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            ),
            previous_names={
                "vpar": "delta_g(mp,var_vpar,ms,mt)",
                "rho": "delta_rho_g",
                "psi": "delta_ps",
            },
        )
    if "rhoimp" in rows:
        equation = impurity_density_equation_rhoimp(with_TiTe=bool(with_tite))
        _add_linearized(
            pools, "rhoimp",
            equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            ),
            previous_names={"rhoimp": "delta_g(mp,var_rhoimp,ms,mt)"},
        )
    if "rhon" in rows:
        equation = neutral_density_equation_rhon(with_TiTe=bool(with_tite))
        _add_linearized(
            pools, "rhon",
            equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            ),
            previous_names={"rhon": "delta_g(mp,var_rhon,ms,mt)"},
        )
    # The ion and electron energy equations exist only in the two-temperature
    # branch of the element routine.
    if "ti" in rows and with_tite:
        equation = ion_energy_equation_Ti(
            st_form=st_form,
        )
        _add_linearized(
            pools, "ti",
            equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            ),
            previous_names={
                "Ti": "delta_g(mp,var_Ti,ms,mt)",
                "rho": "delta_rho_g",
                "rhoimp": "delta_g(mp,var_rhoimp,ms,mt)",
            },
        )
    if "te" in rows and with_tite:
        equation = electron_energy_equation_Te(
            st_form=st_form,
        )
        _add_linearized(
            pools, "te",
            equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            ),
            previous_names={
                "Te": "delta_g(mp,var_Te,ms,mt)",
                "rho": "delta_rho_g",
                "rhoimp": "delta_g(mp,var_rhoimp,ms,mt)",
            },
        )
    # The single-temperature energy equation exists only in the other branch.
    if "t" in rows and not with_tite:
        equation = total_energy_equation_T(st_form=st_form)
        _add_linearized(
            pools, "t",
            equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            ),
            previous_names={
                "T": "delta_g(mp,var_T,ms,mt)",
                "rho": "delta_rho_g",
                "rhoimp": "delta_g(mp,var_rhoimp,ms,mt)",
            },
        )
    if "zj" in rows:
        _add_linearized(
            pools, "zj", current_constraint_equation_zj().linearize(fields=FIELDS)
        )
    if "w" in rows:
        _add_linearized(
            pools, "w", vorticity_constraint_equation_w().linearize(fields=FIELDS)
        )
    if "u" in rows:
        temperature_branches = (False, True) if with_tite is None else (with_tite,)
        for with_tite in temperature_branches:
            equation = momentum_equation_2(
                with_TiTe=with_tite,
                include_neo=include_neo,
            )
            # The source routine contains both the one- and two-temperature
            # RHS branches.  Retain both raw expressions so terms involving
            # ``Ti0``/``Te0`` can be aligned as well as the single-T terms.
            _add_raw_rhs(
                pools, "u", _raw_evolution_rhs(equation, FIELDS, timestep, zeta),
                previous_names={"u": "delta_u", "rho": "delta_rho_g"},
            )
            result = equation.linearize(
                fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
            )
            _add_linearized(
                pools, "u", result,
                previous_names={"u": "delta_u", "rho": "delta_rho_g"},
                include_rhs=False,
            )
    # For the momentum report retain repeated terms: the Fortran routine has
    # separate base/extension assignments whose identical monomials must be
    # available more than once for source-row alignment.  Other rows retain
    # the historical duplicate suppression.
    if rows == {"u"}:
        return {lhs: list(terms) for lhs, terms in pools.items()}
    return {
        lhs: list(dict.fromkeys(terms))
        for lhs, terms in pools.items()
    }





def _expanded_ordered_monomials(text, *, temperature_model=None):
    """Expand generated text using the same temperature convention as source."""

    variants = _parsed_terms(text, temperature_model=temperature_model)
    if not variants:
        return []
    terms = variants[0]
    common_symbols = set.intersection(
        *(set(term.free_symbols) for term in terms)
    ) if terms else set()
    return sorted(
        terms,
        key=lambda term: _source_term_order(text, term, common_symbols),
    )


def _monomial_structure(term):
    """Return the coefficient-free part of a monomial."""

    return sp.cancel(term).as_coeff_Mul()[1]


def _geometric_structure(term):
    """Monomial structure with the cylindrical radius factored out.

    The element routine and the generator sometimes disagree by a power of
    ``BigR`` alone.  Keying on the rest of the monomial lets such a pair be
    shown on one report line instead of drifting apart, which is how a wrong
    power of the radius becomes visible in Meld.
    """

    structure = _monomial_structure(term)
    radius = sp.Symbol("BigR")
    exponent = sp.degree(sp.together(structure).as_numer_denom()[0], radius)
    try:
        return sp.cancel(structure / radius**exponent)
    except Exception:
        return structure


def _scaling_structure(term):
    """Monomial structure with every scalar scaling factor removed.

    On top of the cylindrical radius this drops the implicitness factor
    ``theta``.  Both are pure scalings of an otherwise identical term, so a
    pair that differs only in them belongs on one report line; that is how a
    tangent written without its ``theta`` becomes visible in Meld.
    """

    structure = _geometric_structure(term)
    theta = sp.Symbol("theta")
    exponent = sp.degree(sp.together(structure).as_numer_denom()[0], theta)
    try:
        return sp.cancel(structure / theta**exponent)
    except Exception:
        return structure


def _combine_source_duplicates(source_monomials, generated_by_key):
    """Sum source monomials the generated side reports only once.

    A Fortran assignment writes a Jacobian as a sum of separately factored
    outer terms, and two of those outer terms can expand onto the same
    monomial.  The product rule in the generator, by contrast, always returns
    fully combined monomials, so the same contribution appears there with the
    summed coefficient.  The example in model 600 is the Taylor-Galerkin
    ``tgnum_u`` tangent, where ``0.25*(w0_x*u_y - w0_y*u_x)*(v_x*u0_y - ...)``
    and ``0.25*(w0_x*u0_y - w0_y*u0_x)*(v_x*u_y - ...)`` both contain
    ``v_x*u0_y*u_y*w0_x`` and therefore combine to ``0.5``.

    Source monomials are combined only when the generated block really has
    fewer copies of that monomial structure; otherwise the one-slot-per-outer
    term layout of the report is preserved.
    """

    generated_structures = defaultdict(int)
    for terms in generated_by_key.values():
        for term in terms:
            generated_structures[_monomial_structure(term)] += 1

    positions = defaultdict(list)
    for index, (_, source_term, _) in enumerate(source_monomials):
        if source_term is None:
            continue
        positions[_monomial_structure(source_term)].append(index)

    dropped = set()
    combined = {}
    for structure, indices in positions.items():
        if len(indices) <= max(generated_structures.get(structure, 0), 1):
            continue
        total = sp.Add(*(source_monomials[index][1] for index in indices))
        combined[indices[0]] = sp.expand(total)
        dropped.update(indices[1:])

    result = []
    for index, entry in enumerate(source_monomials):
        if index in dropped:
            continue
        if index in combined:
            if combined[index] == 0:
                continue
            assignment, _, source_raw = entry
            entry = (assignment, combined[index], source_raw)
        result.append(entry)
    return result


def _source_monomial_slots(assignments, generated_terms, *, temperature_model=None):
    """Expand both sides and align every generated monomial to source order."""

    source_monomials = []
    # The NEO psi tangent is written in JOREK as four additive variations,
    # while SymPy combines equal monomials when differentiating the compact
    # residual.  Combine the source block first as well; otherwise a source
    # coefficient such as ``1`` is compared with the algebraically equivalent
    # combined generated coefficient ``3`` on a different source line.
    source_assignments = assignments
    if (
        assignments
        and assignments[0].lhs == "amat(var_u,var_psi)"
        and any("amu_neo_prof" in item.expression for item in assignments)
    ):
        combined = " + ".join(item.expression for item in assignments)
        source_assignments = (
            SourceAssignment(assignments[0].lhs, assignments[0].row,
                             assignments[0].line, combined),
        )
    combine_neo_psi = (
        assignments
        and assignments[0].lhs == "amat(var_u,var_psi)"
        and any("amu_neo_prof" in item.expression for item in assignments)
    )
    for assignment in source_assignments:
        pieces = ([assignment.expression] if combine_neo_psi
                  else _split_top_level(assignment.expression))
        for piece in pieces:
            variants = _parsed_terms(piece, temperature_model=temperature_model)
            use_two_temperature = (
                assignment.lhs in {
                    "amat(var_psi,var_psi)",
                    "amat(var_psi,var_rho)",
                    "amat_n(var_psi,var_rho)",
                }
                or "var_te" in assignment.lhs
            )
            expanded = []
            if variants:
                if use_two_temperature:
                    # Explicit Te0 source terms already identify the branch;
                    # abstract Pe0 terms require the two-temperature alias.
                    variant_index = (
                        0 if "te0" in piece.lower() else len(variants) - 1
                    )
                else:
                    variant_index = 0
                expanded = list(variants[variant_index])
                common_symbols = set.intersection(
                    *(set(term.free_symbols) for term in expanded)
                ) if expanded else set()
                expanded.sort(
                    key=lambda term: _source_term_order(
                        piece, term, common_symbols
                    )
                )
            if not expanded:
                # Keep unsupported source expressions (currently NEO terms)
                # visible in the report instead of silently dropping them.
                source_monomials.append((assignment, None, piece))
                continue
            for term in expanded:
                if term == 0:
                    continue
                source_monomials.append((assignment, term, None))

    generated_by_key = defaultdict(list)
    generated_by_structure = defaultdict(list)
    generated_by_radius = defaultdict(list)
    generated_by_scaling = defaultdict(list)
    # For a source NEO block, keep only generated NEO terms in this local
    # alignment.  The complete generated pool also contains the ordinary
    # momentum terms, which must not appear as unrelated generated-only rows
    # in the NEO subsection.
    neo_block = bool(
        assignments
        and assignments[0].lhs == "amat(var_u,var_psi)"
        and any("amu_neo_prof" in item.expression for item in assignments)
    )
    pool_terms = (
        [text for text in generated_terms if "amu_neo_prof" in text]
        if neo_block else generated_terms
    )
    for text in pool_terms:
        for term in _expanded_ordered_monomials(
            text, temperature_model=temperature_model,
        ):
            # The residual history in the element routine uses the physical
            # background density ``r0``; only the velocity tangent uses the
            # corrected coefficient ``r0_corr``.  The DSL shares one A form,
            # so apply that source convention when matching RHS monomials.
            if assignments and assignments[0].lhs == "rhs_ij(var_u)":
                term = term.xreplace({sp.Symbol("r0_corr"): sp.Symbol("r0")})
            # JOREK stores W_dia_rho as a logarithmic density derivative
            # (the explicit trial-density factor is absorbed in the work
            # variable), whereas the DSL variation exposes that factor.
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_rho)"
                and term.has(sp.Symbol("W_dia_rho"))
                and term.has(sp.Symbol("rho"))
            ):
                term = sp.cancel(term / sp.Symbol("rho"))
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_ti)"
                and term.has(sp.Symbol("W_dia_Ti"))
                and term.has(sp.Symbol("Ti"))
            ):
                term = sp.cancel(term / sp.Symbol("Ti"))
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_t)"
                and term.has(sp.Symbol("W_dia_T"))
                and term.has(sp.Symbol("T"))
            ):
                term = sp.cancel(term / sp.Symbol("T"))
            if (
                assignments
                and assignments[0].lhs == "amat(var_u,var_rhoimp)"
                and term.has(sp.Symbol("alpha_e_bis"))
                and term.has(sp.Symbol("Te0"))
            ):
                # alpha_e_bis is supplied as d(alpha_e*Te)/dTe, so the work
                # variable already contains the explicit background Te0.
                term = sp.cancel(term / sp.Symbol("Te0"))
            generated_by_key[_canonical_monomial(term)].append(term)
            generated_by_structure[_monomial_structure(term)].append(term)
            generated_by_radius[_geometric_structure(term)].append(term)
            generated_by_scaling[_scaling_structure(term)].append(term)

    source_monomials = _combine_source_duplicates(
        source_monomials, generated_by_key
    )

    def _take(pool, key):
        """Pop the first unconsumed term registered under ``key``."""

        while pool[key]:
            term = pool[key].pop(0)
            if id(term) in consumed:
                continue
            consumed.add(id(term))
            return term
        return None

    # Match in three passes rather than term by term.  A single pass lets an
    # early source monomial with no exact partner consume, through one of the
    # relaxed keys, a generated monomial that a later source monomial matches
    # exactly; the pair then drifts apart in the report even though the block
    # agrees.  Exact matches are therefore all resolved first, then the ones
    # that differ only by a coefficient, then the ones that differ by a power
    # of the cylindrical radius as well.
    consumed = set()
    entries = []
    for assignment, source_term, source_raw in source_monomials:
        if source_term is None:
            entries.append([assignment, source_raw, "", None])
            continue
        entries.append([
            assignment, _canonical_monomial(source_term), "", source_term,
        ])
    for pool, key_of in (
        (generated_by_key, _canonical_monomial),
        (generated_by_structure, _monomial_structure),
        (generated_by_radius, _geometric_structure),
        (generated_by_scaling, _scaling_structure),
    ):
        for entry in entries:
            source_term = entry[3]
            if source_term is None or entry[2]:
                continue
            match = _take(pool, key_of(source_term))
            if match is not None:
                entry[2] = _canonical_monomial(match)
    slots = [(assignment, source, generated)
             for assignment, source, generated, _ in entries]

    # Do not hide generated-only terms.  They are emitted after the aligned
    # source rows with an empty source column, which makes sign or convention
    # inconsistencies visible in Meld.
    synthetic = SourceAssignment(
        assignments[0].lhs, assignments[0].row, 0, ""
    ) if assignments else None
    if synthetic is not None:
        emitted = set()
        for values in generated_by_key.values():
            for term in values:
                if id(term) in consumed:
                    continue
                rendered = _canonical_monomial(term)
                if rendered in emitted:
                    continue
                emitted.add(rendered)
                slots.append((synthetic, "", rendered))
    return slots


def _unmatched_generated_blocks(rows, source_lhs, generated, temperature_model):
    """Emit generated assignments the source routine does not write at all.

    A Jacobian column that is missing from the element routine entirely has no
    source assignment to align against, so it would otherwise disappear from
    the report instead of showing up as an unimplemented block.
    """

    slots = []
    for lhs, terms in generated.items():
        if lhs in source_lhs or not terms:
            continue
        row_match = re.search(r"\(var_(psi|u|zj|w|rho|rhoimp|vpar|ti|te|t)\b", lhs)
        if row_match is None or row_match.group(1) not in rows:
            continue
        synthetic = SourceAssignment(lhs, row_match.group(1), 0, "")
        emitted = set()
        for text in terms:
            for term in _expanded_ordered_monomials(
                text, temperature_model=temperature_model,
            ):
                rendered = _canonical_monomial(term)
                if rendered in emitted:
                    continue
                emitted.add(rendered)
                slots.append((synthetic, "", rendered))
    return slots


def _source_term_order(source_piece, term, common_symbols=()):
    """Estimate textual order for expanded terms inside one source term."""

    text = source_piece.lower()
    positions = []
    common_symbols = set(common_symbols)
    for symbol in term.free_symbols:
        if symbol in common_symbols:
            continue
        name = str(symbol).lower()
        # A symbol appearing once is usually the factor that distinguishes this
        # expanded monomial (zj0/current_source/Jb in the induction example).
        if len(re.findall(r"\b{}\b".format(re.escape(name)), text)) == 1:
            positions.append(text.find(name))
    return min(positions) if positions else len(text)


def _factor_order(factor):
    """Use a stable JOREK-oriented order for multiplicative factors."""

    base = factor
    while getattr(base, "is_Pow", False):
        base = base.base
    name = str(base)
    prefixes = (
        "F0", "BigR", "Z", "v", "u", "psi", "zj", "w", "rho",
        "T", "Ti", "Te", "xjac", "theta", "tstep", "zeta",
    )
    prefix_rank = next(
        (index for index, prefix in enumerate(prefixes) if name.startswith(prefix)),
        len(prefixes),
    )
    return prefix_rank, name, str(factor)


def _canonical_monomial(term):
    """Render one monomial with the same factor order in both reports."""

    # SymPy may leave a common numeric factor in both the numerator and an
    # expanded NEO denominator (for example ``-2/(2*D)``).  Reduce rational
    # factors before ordering products so algebraically identical source and
    # generated monomials receive the same key.
    term = sp.cancel(term)
    coefficient = sp.S.One
    factors = []
    for factor in sp.Mul.make_args(term):
        if factor.is_Number:
            coefficient *= factor
        else:
            factors.append(factor)
    factors.sort(key=_factor_order)
    numerator = [factor for factor in factors if not (factor.is_Pow and factor.exp.is_negative)]
    denominator = [factor.base ** (-factor.exp) for factor in factors if factor.is_Pow and factor.exp.is_negative]
    pieces = []
    if coefficient not in (1, -1) or (not numerator and not denominator):
        pieces.append(sp.sstr(coefficient))
    pieces.extend(sp.sstr(factor) for factor in numerator)
    result = "*".join(pieces) if pieces else "1"
    for factor in denominator:
        rendered = sp.sstr(factor)
        if not factor.is_Atom:
            rendered = "({})".format(rendered)
        result += "/{}".format(rendered)
    if coefficient == -1 and (numerator or denominator):
        result = "-{}".format(result)
    return result


def _canonical_display(text, *, temperature_model=None):
    """Normalize multiplication ordering without changing term grouping.

    The temperature convention must be named rather than picked by position:
    ``_parsed_terms`` drops variants that coincide, so an expression without
    two-temperature aliases would otherwise be displayed through whichever
    convention happened to survive.
    """

    variants = _parsed_terms(text, temperature_model=temperature_model)
    if not variants:
        return text
    terms = variants[0]
    if len(terms) == 1 and terms[0] == 0:
        return ""
    return " + ".join(_canonical_monomial(term) for term in terms)


def _slot_mismatches(slots):
    """Count report lines that are blank on one side."""

    return sum(1 for _, source, generated in slots if not source or not generated)


def _best_slots(lhs_assignments, pools, temperature_model):
    """Align against whichever generated spelling fits the source block.

    The element routine writes a poloidal bracket sometimes as
    ``a_s*b_t - a_t*b_s`` and sometimes as ``xjac*(a_x*b_y - a_y*b_x)``, and
    it picks differently for a residual and for one of its own tangents.  The
    two spellings are identical, but they expand into different monomials, so
    aligning them line by line requires using the spelling the block at hand
    happens to use.
    """

    best = None
    for pool in pools:
        slots = _source_monomial_slots(
            lhs_assignments, pool, temperature_model=temperature_model,
        )
        score = _slot_mismatches(slots)
        if best is None or score < best[0]:
            best = (score, slots)
        if score == 0:
            break
    return best[1]


def _source_slots(assignments, generated, *, temperature_model=None,
                  generated_alternative=None):
    """Align generated monomials to source monomials in source order."""

    alternative = generated_alternative or {}

    slots = []
    used = defaultdict(set)
    positional_consumed = defaultdict(int)
    # For a model-600 momentum report, use the same monomial alignment for
    # every RHS/AMAT block.  Each assignment is matched independently because
    # the source contains separate columns and toroidal channels.  This keeps
    # the generated line in the same row as its source monomial and leaves a
    # genuinely empty line when that contribution is not implemented.
    if assignments and all(assignment.row == "u" for assignment in assignments):
        def _remove_branch_duplicates(pool):
            """Drop exact terms repeated only because both T branches run.

            Terms containing an explicit background temperature belong to a
            branch and must remain separate.  Terms without ``T0``, ``Ti0``
            or ``Te0`` are branch-independent (base, magnetic, and most
            geometric terms); retaining two identical copies made the report
            look as if the generator had duplicated physics.
            """
            seen = set()
            result = []
            for term in pool:
                lowered = term.lower()
                branch_specific = any(
                    marker in lowered for marker in ("t0", "ti0", "te0")
                )
                if not branch_specific:
                    if term in seen:
                        continue
                    seen.add(term)
                result.append(term)
            return result

        by_lhs = defaultdict(list)
        for assignment in assignments:
            by_lhs[assignment.lhs].append(assignment)
        slots = []
        for lhs, lhs_assignments in by_lhs.items():
            pool = generated.get(lhs, ())
            neo_assignments = [
                item for item in lhs_assignments
                if "amu_neo_prof" in item.expression
            ]
            regular_assignments = [
                item for item in lhs_assignments
                if "amu_neo_prof" not in item.expression
            ]
            if regular_assignments:
                slots.extend(
                    _source_monomial_slots(
                        regular_assignments,
                        _remove_branch_duplicates(
                            [term for term in pool if "amu_neo_prof" not in term]
                        ), temperature_model=temperature_model,
                    )
                )
            if neo_assignments:
                slots.extend(
                    _source_monomial_slots(
                        neo_assignments,
                        [term for term in pool if "amu_neo_prof" in term],
                        temperature_model=temperature_model,
                    )
                )
        slots.extend(_unmatched_generated_blocks(
            {"u"}, set(by_lhs), generated, temperature_model,
        ))
        return slots
    # Export the remaining equations (notably PSI) at monomial granularity as
    # well.  The older path below aligned whole Fortran additive pieces, which
    # made the PSI report look as though monomials were missing even when the
    # generated expression contained them.
    if assignments:
        by_lhs = defaultdict(list)
        for assignment in assignments:
            by_lhs[assignment.lhs].append(assignment)
        slots = [
            slot
            for lhs, lhs_assignments in by_lhs.items()
            for slot in _best_slots(
                lhs_assignments,
                [generated.get(lhs, ())]
                + ([alternative[lhs]] if lhs in alternative else []),
                temperature_model,
            )
        ]
        slots.extend(_unmatched_generated_blocks(
            {assignment.row for assignment in assignments},
            set(by_lhs), generated, temperature_model,
        ))
        return slots
    # ``assignments`` is empty only when the requested row has no source
    # assignment in this temperature branch; the generator produces nothing
    # for it either, so there is nothing to align.
    return []


def _render_report(
    slots, side, rows=ROWS, *, include_neo=False, temperature_model=None,
):
    section = {
        "single": "Single-temperature (T) model",
        "two": "Two-temperature (Ti/Te) model",
    }.get(temperature_model)
    lines = [
        "## {}".format(section) if section else "# Model 600 equation terms",
        "",
        "> Rows currently exported: {}.".format(
            ", ".join("`{}`".format(row) for row in rows)
        ),
        (
            "> NEO RHS and AMAT terms are included when the corresponding source branch is enabled."
            if include_neo
            else "> NEO RHS and AMAT terms are omitted from this comparison."
        ),
        "",
    ]
    for row in rows:
        lines.extend(("### var_{}".format(row) if section else "## var_{}".format(row), ""))
        row_slots = [slot for slot in slots if slot[0].row == row]
        lhs_order = list(dict.fromkeys(slot[0].lhs for slot in row_slots))
        for lhs in lhs_order:
            lines.extend(("#### `{}`".format(lhs) if section else "### `{}`".format(lhs), ""))
            for assignment, source, generated in row_slots:
                if assignment.lhs != lhs:
                    continue
                value = source if side == "source" else generated
                value = " ".join(value.split())
                value = _canonical_display(
                    value, temperature_model=temperature_model,
                )
                # Keep exactly one physical line per aligned term.  An
                # unavailable term is intentionally an actually empty line,
                # which makes its position immediately visible in Meld.
                lines.append(value)
            lines.append("")
    return "\n".join(lines) + "\n"


def export_model600_markdown(
    source_path, source_output, generated_output, equations=None,
    include_neo=False,
):
    """Write aligned source/generated reports and return their paths.

    ``equations`` may contain any of ``psi``, ``u``, ``zj`` and ``w``.  If it
    is omitted, all four rows are exported.  NEO terms are omitted by default;
    pass ``include_neo=True`` to include them.
    """

    # The Fortran spells the temperature rows ``var_Ti``/``var_Te``; the
    # exporter keys everything on the lower-case assignment name.
    rows = tuple(row.lower() for row in (equations or ROWS))
    invalid = set(rows) - set(ROWS)
    if invalid:
        raise ValueError("Unknown model-600 equation(s): {}".format(", ".join(sorted(invalid))))
    all_assignments = _read_source_assignments(source_path)
    source_sections = []
    generated_sections = []
    for temperature_model, with_tite in (("single", False), ("two", True)):
        source_assignments = [
            assignment for assignment in all_assignments
            if assignment.row in rows
            and assignment.temperature_model in ("both", temperature_model)
            and (include_neo or "amu_neo_prof" not in assignment.expression)
        ]
        generated = _generated_pools(
            rows, include_neo=include_neo,
            with_tite=with_tite,
        )
        generated_alternative = {}
        if {"vpar", "ti", "te", "t"} & set(rows):
            # The parallel-velocity row is the one whose source assignments
            # disagree among themselves about how to spell a poloidal
            # bracket, so build the other spelling as well.
            alternative_rows = tuple(
                row for row in ("vpar", "ti", "te", "t") if row in rows
            )
            alternative_pools = _generated_pools(
                alternative_rows,
                include_neo=include_neo, with_tite=with_tite,
                st_form=False,
            )
            generated_alternative = dict(alternative_pools)
        slots = _source_slots(
            source_assignments, generated, temperature_model=temperature_model,
            generated_alternative=generated_alternative,
        )
        source_sections.append(
            _render_report(
                slots, "source", rows, include_neo=include_neo,
                temperature_model=temperature_model,
            )
        )
        generated_sections.append(
            _render_report(
                slots, "generated", rows, include_neo=include_neo,
                temperature_model=temperature_model,
            )
        )
    source_output = Path(source_output)
    generated_output = Path(generated_output)
    source_output.parent.mkdir(parents=True, exist_ok=True)
    generated_output.parent.mkdir(parents=True, exist_ok=True)
    header = "# Model 600 equation terms\n\n"
    source_output.write_text(header + "\n".join(source_sections), encoding="utf-8")
    generated_output.write_text(header + "\n".join(generated_sections), encoding="utf-8")
    return source_output, generated_output
