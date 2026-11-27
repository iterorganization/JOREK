"""Reading the JOREK element routine and comparing its terms with the DSL.

This module owns step 1 and step 3 of the pipeline: it extracts the assembled
``rhs_ij``/``amat`` expressions from a model's ``mod_elt_matrix_fft.f90``,
rewrites the element routine's work variables into the names the symbolic DSL
prints, parses the result into SymPy, and compares two such maps.
"""

from pathlib import Path
import re
from collections import Counter
from typing import Dict, Iterable, Mapping

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

from .channels import split_toroidal_channels
from .fortran import fortran
from .model199 import (
    FIELDS,
    T,
    current_constraint_equation_3,
    density_equation_5,
    induction_equation_1,
    j,
    momentum_equation_2,
    psi,
    rho,
    u,
    vorticity_constraint_equation_4,
    temperature_equation_6,
    xjac,
)
from .operators import expand_brackets, lower_brackets_to_element
from .symbols import coefficient


EQUATION_1_ASSIGNMENTS = (
    "rhs_ij_1",
    "amat_11",
    "amat_12",
    "amat_12_n",
    "amat_13",
    "amat_16",
)

EQUATION_2_ASSIGNMENTS = (
    "rhs_ij_2", "amat_21", "amat_22", "amat_23", "amat_23_n",
    "amat_24", "amat_25", "amat_26",
)
EQUATION_3_ASSIGNMENTS = ("rhs_ij_3", "amat_31", "amat_33")
EQUATION_4_ASSIGNMENTS = ("rhs_ij_4", "amat_42", "amat_44")
EQUATION_5_ASSIGNMENTS = (
    "rhs_ij_5", "rhs_ij_5_k", "amat_51", "amat_51_k", "amat_52",
    "amat_55", "amat_55_k", "amat_55_n", "amat_55_kn",
)
EQUATION_6_ASSIGNMENTS = (
    "rhs_ij_6", "rhs_ij_6_k", "amat_61", "amat_61_k", "amat_62",
    "amat_63", "amat_66", "amat_66_k", "amat_66_n", "amat_66_kn",
)


class FortranComparisonError(AssertionError):
    """Generated and source Fortran expressions differ."""


def format_expression_terms(expression) -> str:
    """Format an expression as one additive term per line.

    The comparison itself remains symbolic; this is only a diagnostic
    representation.  Expanding the outer addition makes it much easier to
    spot a missing sign or coefficient in the generated weak-form terms.
    """

    expanded = sp.expand(expression)
    if expanded == 0:
        return "    0"
    terms = sp.Add.make_args(expanded)
    lines = []
    for term in terms:
        if term.could_extract_minus_sign():
            lines.append("    - {}".format(str(-term)))
        else:
            lines.append("    + {}".format(str(term)))
    return "\n".join(lines)


def format_assignment(name: str, expression, label: str = "") -> str:
    """Return a readable, multiline assignment diagnostic."""

    heading = "{}{}:".format(label + " " if label else "", name)
    return heading + "\n" + format_expression_terms(expression)


def format_term_difference(source, generated) -> str:
    """Show unmatched additive terms without subtracting the expressions."""

    source_terms = Counter(sp.Add.make_args(sp.expand(source)))
    generated_terms = Counter(sp.Add.make_args(sp.expand(generated)))
    source_only = list((source_terms - generated_terms).elements())
    generated_only = list((generated_terms - source_terms).elements())
    lines = []
    if source_only:
        lines.append("    source-only:")
        lines.extend("      {}".format(_term_text(term)) for term in source_only)
    if generated_only:
        lines.append("    generated-only:")
        lines.extend("      {}".format(_term_text(term)) for term in generated_only)
    return "\n".join(lines) if lines else "    (no unmatched terms)"


def _term_text(term) -> str:
    """Render one term while retaining its sign."""

    if term.could_extract_minus_sign():
        return "- {}".format(str(-term))
    return "+ {}".format(str(term))


def extract_fortran_assignments(path, names: Iterable[str]) -> Dict[str, str]:
    """Extract continued scalar assignments from a Fortran source file."""

    requested = tuple(names)
    patterns = {
        name: re.compile(r"^\s*{}\s*=\s*(.*)$".format(re.escape(name)), re.I)
        for name in requested
    }
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    found = {}
    index = 0
    while index < len(lines):
        uncommented = lines[index].split("!", 1)[0]
        matched_name = None
        matched_rhs = None
        for name, pattern in patterns.items():
            match = pattern.match(uncommented)
            if match:
                matched_name = name
                matched_rhs = match.group(1)
                break
        if matched_name is None:
            index += 1
            continue

        pieces = []
        current = matched_rhs
        while True:
            current = current.strip()
            continued = current.endswith("&")
            pieces.append(current[:-1] if continued else current)
            if not continued:
                break
            index += 1
            if index >= len(lines):
                raise ValueError("Unterminated assignment for {}".format(matched_name))
            current = lines[index].split("!", 1)[0].lstrip()
            if current.startswith("&"):
                current = current[1:]
        found[matched_name] = " ".join(piece.strip() for piece in pieces)
        index += 1

    missing = [name for name in requested if name not in found]
    if missing:
        raise ValueError(
            "Assignments not found in {}: {}".format(path, ", ".join(missing))
        )
    return found


_D_LITERAL = re.compile(r"(?i)(?<![A-Za-z_])(\d+(?:\.\d*)?|\.\d+)d([+-]?\d+)")
_IDENTIFIER = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")


def normalize_fortran_text(expression: str) -> str:
    """Normalize the small Fortran-expression subset used by the comparator."""

    result = expression.replace("&", " ")
    result = re.sub(r"\bBigr\b", "BigR", result, flags=re.I)
    result = re.sub(
        r"current_source\s*\(\s*ms\s*,\s*mt\s*\)",
        "current_source",
        result,
        flags=re.I,
    )
    for function_name in ("particle_source", "heat_source_i", "heat_source_e",
                          "heat_source"):
        result = re.sub(
            r"{}\s*\(\s*ms\s*,\s*mt\s*\)".format(function_name),
            function_name,
            result,
            flags=re.I,
        )
    result = re.sub(
        r"delta_g\s*\(\s*mp\s*,\s*1\s*,\s*ms\s*,\s*mt\s*\)",
        "delta_g_1",
        result,
        flags=re.I,
    )
    result = re.sub(
        r"delta_g\s*\(\s*mp\s*,\s*([56])\s*,\s*ms\s*,\s*mt\s*\)",
        r"delta_g_\1",
        result,
        flags=re.I,
    )
    aliases = {
        "D_par_local": "D_par",
        "vv2": "(BigR**2*(u0_x**2+u0_y**2))",
        "r0_hat": "(BigR**2*r0)",
        "r0_x_hat": "(2*BigR*BigR_x*r0+BigR**2*r0_x)",
        "r0_y_hat": "(BigR**2*r0_y)",
        "rho_hat": "(BigR**2*rho)",
        "rho_x_hat": "(2*BigR*BigR_x*rho+BigR**2*rho_x)",
        "rho_y_hat": "(BigR**2*rho_y)",
        "P0_s": "(r0_s*T0+r0*T0_s)",
        "P0_t": "(r0_t*T0+r0*T0_t)",
        "BB2": "((F0**2+ps0_x**2+ps0_y**2)/BigR**2)",
        "BB2_psi": "(2*(psi_x*ps0_x+psi_y*ps0_y)/BigR**2)",
        "Bgrad_rho_star": "((v_x*ps0_y-v_y*ps0_x)/BigR)",
        "Bgrad_rho_k_star": "(F0*v_p/BigR**2)",
        "Bgrad_rho": "((F0*r0_p/BigR+r0_x*ps0_y-r0_y*ps0_x)/BigR)",
        "Bgrad_rho_star_psi": "((v_x*psi_y-v_y*psi_x)/BigR)",
        "Bgrad_rho_psi": "((r0_x*psi_y-r0_y*psi_x)/BigR)",
        "Bgrad_rho_rho": "((rho_x*ps0_y-rho_y*ps0_x)/BigR)",
        "Bgrad_rho_rho_n": "(F0*rho_p/BigR**2)",
        "Bgrad_T_star": "((v_x*ps0_y-v_y*ps0_x)/BigR)",
        "Bgrad_T_k_star": "(F0*v_p/BigR**2)",
        "Bgrad_T": "((F0*T0_p/BigR+T0_x*ps0_y-T0_y*ps0_x)/BigR)",
        "Bgrad_T_star_psi": "((v_x*psi_y-v_y*psi_x)/BigR)",
        "Bgrad_T_psi": "((T0_x*psi_y-T0_y*psi_x)/BigR)",
        "Bgrad_T_T": "((T_x*ps0_y-T_y*ps0_x)/BigR)",
        "Bgrad_T_T_n": "(F0*T_p/BigR**2)",
    }
    for name, replacement in aliases.items():
        result = re.sub(r"\b{}\b".format(name), replacement, result, flags=re.I)
    result = _D_LITERAL.sub(lambda match: match.group(1) + "e" + match.group(2), result)
    return " ".join(result.split())


def parse_fortran_expression(expression: str):
    """Parse a normalized scalar Fortran expression into an exact SymPy tree."""

    normalized = normalize_fortran_text(expression)
    names = set(_IDENTIFIER.findall(normalized))
    local_dict = {name: sp.Symbol(name) for name in names}
    parsed = parse_expr(normalized, local_dict=local_dict, evaluate=True)
    return sp.nsimplify(parsed)


def generated_model199_equation1_assignments() -> Dict[str, sp.Expr]:
    """Generate equation-1 terms in the same local convention as the source."""

    equation = induction_equation_1()
    result = equation.linearize(
        fields=FIELDS,
        timestep=coefficient("tstep"),
        theta=coefficient("theta"),
        zeta=coefficient("zeta"),
    )

    rhs_channels = split_toroidal_channels(result.rhs)
    psi_channels = split_toroidal_channels(result.amat[psi])
    u_channels = split_toroidal_channels(result.amat[u])
    j_channels = split_toroidal_channels(result.amat[j])
    T_channels = split_toroidal_channels(result.amat[T])

    unsupported = {
        "rhs": (rhs_channels.n, rhs_channels.k, rhs_channels.kn),
        "psi": (psi_channels.n, psi_channels.k, psi_channels.kn),
        "u": (u_channels.k, u_channels.kn),
        "j": (j_channels.n, j_channels.k, j_channels.kn),
        "T": (T_channels.n, T_channels.k, T_channels.kn),
    }
    for label, values in unsupported.items():
        if any(value != 0 for value in values):
            raise ValueError(
                "Unexpected toroidal channel in equation-1 {} block".format(label)
            )

    generated = {
        "rhs_ij_1": rhs_channels.p,
        "amat_11": psi_channels.p,
        "amat_12": u_channels.p,
        "amat_12_n": u_channels.n,
        "amat_13": j_channels.p,
        "amat_16": T_channels.p,
    }
    return {
        name: lower_brackets_to_element(expression, xjac)
        for name, expression in generated.items()
    }


def generated_assignment_text() -> Dict[str, str]:
    """Render generated terms with the identifiers used in model 199."""

    return {
        name: fortran(
            expression,
            previous_names={"psi": "delta_g(mp,1,ms,mt)"},
        )
        for name, expression in generated_model199_equation1_assignments().items()
    }


def generated_partial_model199_assignments() -> Dict[str, str]:
    """Render generated assignments for the currently implemented equations."""

    timestep = coefficient("tstep")
    theta = coefficient("theta")
    zeta = coefficient("zeta")
    linearized = momentum_equation_2().linearize(
        fields=FIELDS, timestep=timestep, theta=theta, zeta=zeta
    )
    channels = split_toroidal_channels(linearized.rhs)
    result_map = {"rhs_ij_2": channels.p}
    for name, field, channel in (
        ("amat_21", psi, "p"), ("amat_22", u, "p"),
        ("amat_23", j, "p"), ("amat_23_n", j, "n"),
        ("amat_24", FIELDS[3], "p"), ("amat_25", FIELDS[4], "p"),
        ("amat_26", T, "p"),
    ):
        result_map[name] = getattr(split_toroidal_channels(linearized.amat[field]), channel)
    for equation, rhs_name, block_fields in (
        (current_constraint_equation_3(), "rhs_ij_3", (("amat_31", psi), ("amat_33", j))),
        (vorticity_constraint_equation_4(), "rhs_ij_4", (("amat_42", u), ("amat_44", FIELDS[3]))),
    ):
        linearized = equation.linearize(fields=FIELDS)
        result_map[rhs_name] = linearized.rhs
        for assignment, field in block_fields:
            result_map[assignment] = linearized.amat[field]
    return {
        name: fortran(value, previous_names={"u": "delta_u"})
        for name, value in result_map.items()
    }


def generated_model199_transport_assignments() -> Dict[str, str]:
    """Render generated RHS/AMAT assignments for equations 5 and 6."""

    timestep = coefficient("tstep")
    theta = coefficient("theta")
    zeta = coefficient("zeta")
    output = {}
    for equation, fields, names in (
        (density_equation_5(), FIELDS, EQUATION_5_ASSIGNMENTS),
        (temperature_equation_6(), FIELDS, EQUATION_6_ASSIGNMENTS),
    ):
        linearized = equation.linearize(
            fields=fields, timestep=timestep, theta=theta, zeta=zeta
        )
        rhs = split_toroidal_channels(linearized.rhs)
        output[names[0]] = rhs.p
        output[names[1]] = rhs.k
        row_blocks = (("51", psi), ("52", u), ("55", rho)) if equation.name.endswith("density") else (("61", psi), ("62", u), ("63", j), ("66", T))
        for block, field in row_blocks:
            blocks = split_toroidal_channels(linearized.amat[field])
            output["amat_{}".format(block)] = blocks.p
            for channel in ("k", "n", "kn"):
                if blocks.__getattribute__(channel) != 0:
                    output["amat_{}_{}".format(block, channel)] = getattr(blocks, channel)
    previous_names = {
        "rho": "delta_g(mp,5,ms,mt)",
        "T": "delta_g(mp,6,ms,mt)",
    }
    return {
        name: fortran(expand_brackets(value), previous_names=previous_names)
        for name, value in output.items()
    }


def compare_assignment_maps(
    source: Mapping[str, str],
    generated: Mapping[str, str],
    *,
    set_BigR_x_one: bool = True,
) -> None:
    """Compare assignment maps and raise with symbolic differences.

    ``set_BigR_x_one`` reflects the usual cylindrical-coordinate convention
    and removes element-geometry ``BigR_x`` factors from this comparison.
    """

    messages = []
    for name in source:
        if name not in generated:
            messages.append("{}: missing generated assignment".format(name))
            continue
        source_expr = parse_fortran_expression(source[name])
        generated_expr = parse_fortran_expression(generated[name])
        if set_BigR_x_one:
            bigr_x = sp.Symbol("BigR_x")
            source_expr = source_expr.subs(bigr_x, 1)
            generated_expr = generated_expr.subs(bigr_x, 1)
        difference = sp.simplify(generated_expr - source_expr)
        if difference != 0:
            messages.append(
                "{name}:\n"
                "{difference}".format(
                    name=name,
                    difference=format_term_difference(source_expr, generated_expr),
                )
            )
    extra = sorted(set(generated) - set(source))
    if extra:
        messages.append("unexpected generated assignments: {}".format(", ".join(extra)))
    if messages:
        raise FortranComparisonError(
            "Equation comparison failed:\n" + "\n".join(messages)
        )


def compare_model199_equation1(path, *, set_BigR_x_one: bool = True) -> None:
    """Compare generated equation 1 with model199/mod_elt_matrix_fft.f90."""

    source = extract_fortran_assignments(path, EQUATION_1_ASSIGNMENTS)
    compare_assignment_maps(
        source, generated_assignment_text(), set_BigR_x_one=set_BigR_x_one
    )


def normalize_model600_text(
    expression, *, single_temperature=False, two_temperature=False
):
    # The element routine uses short background-density names, whereas the
    # symbolic model uses the descriptive field names.  They denote the same
    # frozen quantities in the single-temperature momentum equation.
    expression = re.sub(r"\brn0(?=\b|_)", "rhon0", expression, flags=re.I)
    expression = re.sub(r"\brimp0(?=\b|_)", "rhoimp0", expression, flags=re.I)
    # Convert the physical-coordinate cross product used by the source to its
    # equivalent element-coordinate bracket before expanding aliases.
    expression = re.sub(
        r"ps0_x\s*\*\s*Pe0_y\s*-\s*ps0_y\s*\*\s*Pe0_x",
        "(ps0_s*Pe0_t-ps0_t*Pe0_s)/xjac",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"factor\(\s*var_psi\s*,\s*([1-6])\s*\)",
        "1",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"factor\(\s*var_(?:zj|w)\s*,\s*([1-6])\s*\)",
        "1",
        expression,
        flags=re.I,
    )
    # ``delta_g(mp,var_X,ms,mt)`` is the previous Newton increment of field X.
    # Give every field the same short spelling the DSL prints.
    expression = re.sub(
        r"delta_g\s*\(\s*mp\s*,\s*var_([a-z0-9_]+)\s*,\s*ms\s*,\s*mt\s*\)",
        r"delta_\1",
        expression,
        flags=re.I,
    )
    expression = re.sub(
        r"current_source\(\s*ms\s*,\s*mt\s*\)",
        "current_source",
        expression,
        flags=re.I,
    )
    # NEO profile work arrays are scalar values at the element/toroidal
    # point, despite being called as indexed Fortran functions.
    for profile in ("amu_neo_prof", "aki_neo_prof"):
        expression = re.sub(
            r"\b{}\s*\(\s*ms\s*,\s*mt\s*\)".format(profile),
            profile,
            expression,
            flags=re.I,
        )
    # The symbolic DSL prints the supplied density correction as a function
    # call, while the element routine stores its value and derivative in the
    # work variables ``r0_corr`` and ``dr0_corr_dn``.
    expression = re.sub(
        r"corr_neg_dens\(\s*r0\s*\)", "r0_corr", expression, flags=re.I
    )
    expression = re.sub(r"\bcorr_neg_dens\b", "r0_corr", expression, flags=re.I)
    expression = re.sub(
        r"dr0_corr_dn\(\s*r0\s*\)", "dr0_corr_dn", expression, flags=re.I
    )
    # The element routine's displayed tangent is written for the default
    # correction branch where d(r0_corr)/d(rho)=1; compare that convention
    # with the explicit DSL derivative call.
    expression = re.sub(r"\bdr0_corr_dn\b", "1", expression)
    if single_temperature:
        # The one-temperature model evolves the total temperature.  Do not
        # identify either Ti0 or Te0 with T0 individually; only their sum
        # (and the corresponding directional/spatial derivatives) is T0.
        for suffix in ("", "_s", "_t", "_p", "_x", "_y", "_xx", "_yy", "_xy"):
            expression = re.sub(
                r"\bTi0{}\s*\+\s*Te0{}\b".format(suffix, suffix),
                "T0{}".format(suffix), expression, flags=re.I,
            )
            expression = re.sub(
                r"\bTe0{}\s*\+\s*Ti0{}\b".format(suffix, suffix),
                "T0{}".format(suffix), expression, flags=re.I,
            )
        expression = re.sub(
            r"\balpha_imp_T\b", "(alpha_imp*T0)", expression, flags=re.I,
        )
    # ``alpha_e_T`` is the DSL print name of the closure value alpha_e*Te0.
    # It is spelled the same way in both temperature models; the single
    # temperature substitution Te0 = T0/2 is applied at the end.
    expression = re.sub(
        r"\balpha_e_T\b", "(alpha_e*Te0)", expression, flags=re.I,
    )
    pressure_temperature = "Te0"
    ion_temperature = "Ti0"
    # ``construct_pressure`` always includes impurity pressure when the
    # impurity extension is active.  Reports compare that full definition,
    # rather than treating P0/Pi0/Pe0 as opaque element work variables.
    # ``construct_pressure`` is written once for both temperature models: it
    # always builds the species pressures from Ti0/Te0 and the per-species
    # impurity coefficients.  The one-temperature branch reaches it with
    # Ti0=Te0=T0/2, which the substitutions at the end of this function
    # apply.  Expanding the aliases in the two-species basis therefore
    # reproduces the element routine exactly in both branches.
    ion_alpha = "alpha_i"
    ion_density = "(r0+rhoimp0*alpha_i)"
    electron_density = "(r0+rhoimp0*alpha_e)"
    electron_density_tangent = "(r0+rhoimp0*alpha_e_bis)"

    pi0 = "({}*{})".format(ion_density, ion_temperature)
    pe0 = "({}*{})".format(electron_density, pressure_temperature)
    pi_derivatives = {
        "s": "((r0_s+rhoimp0_s*{a})*{t}+{d}*{t}_s)",
        "t": "((r0_t+rhoimp0_t*{a})*{t}+{d}*{t}_t)",
        "p": "((r0_p+rhoimp0_p*{a})*{t}+{d}*{t}_p)",
        "x": "((r0_x+rhoimp0_x*{a})*{t}+{d}*{t}_x)",
        "y": "((r0_y+rhoimp0_y*{a})*{t}+{d}*{t}_y)",
        "xx": "((r0_xx+rhoimp0_xx*{a})*{t}+2*(r0_x+rhoimp0_x*{a})*{t}_x+{d}*{t}_xx)",
        "yy": "((r0_yy+rhoimp0_yy*{a})*{t}+2*(r0_y+rhoimp0_y*{a})*{t}_y+{d}*{t}_yy)",
        "xy": "((r0_xy+rhoimp0_xy*{a})*{t}+(r0_x+rhoimp0_x*{a})*{t}_y+(r0_y+rhoimp0_y*{a})*{t}_x+{d}*{t}_xy)",
    }
    pi_derivatives = {
        key: value.format(a=ion_alpha, d=ion_density, t=ion_temperature)
        for key, value in pi_derivatives.items()
    }
    pe_derivatives = {
        "s": "((r0_s+rhoimp0_s*alpha_e)*{t}+{dt}*{t}_s)",
        "t": "((r0_t+rhoimp0_t*alpha_e)*{t}+{dt}*{t}_t)",
        "p": "((r0_p+rhoimp0_p*alpha_e)*{t}+{dt}*{t}_p)",
        "x": "((r0_x+rhoimp0_x*alpha_e)*{t}+{dt}*{t}_x)",
        "y": "((r0_y+rhoimp0_y*alpha_e)*{t}+{dt}*{t}_y)",
    }
    pe_derivatives = {
        key: value.format(d=electron_density, dt=electron_density_tangent,
                          t=pressure_temperature)
        for key, value in pe_derivatives.items()
    }
    p0 = "({}+{})".format(pi0, pe0)
    p0_derivatives = {
        key: "({}+{})".format(pi_derivatives[key], pe_derivatives[key])
        for key in ("s", "t", "p")
    }
    # The implicit heating floor of the ion and electron energy equations is
    # written inline with ``min`` and ``exp``.  Name the two pieces so that
    # the expression parser sees ordinary work values, and resolve the
    # derivative of the exponential through its own value.
    expression = re.sub(
        r"\bmin\s*\(\s*Ti0\s*,\s*Tie_min_neg\s*\)", "Ti0_floor",
        expression, flags=re.I,
    )
    expression = re.sub(
        r"\bexp\s*\(\s*\(\s*Ti0_floor\s*-\s*Tie_min_neg\s*\)\s*/\s*"
        r"\(\s*0\.5(?:d0)?\s*\*\s*Tie_min_neg\s*\)\s*\)",
        "Ti_floor_exp", expression, flags=re.I,
    )
    expression = re.sub(
        r"\bdTi_floor_exp\b", "(Ti_floor_exp/(0.5*Tie_min_neg))",
        expression, flags=re.I,
    )
    expression = re.sub(r"\bdTi_floor\b", "1", expression, flags=re.I)
    expression = re.sub(
        r"\bmin\s*\(\s*Te0\s*,\s*Tie_min_neg\s*\)", "Te0_floor",
        expression, flags=re.I,
    )
    expression = re.sub(
        r"\bexp\s*\(\s*\(\s*Te0_floor\s*-\s*Tie_min_neg\s*\)\s*/\s*"
        r"\(\s*0\.5(?:d0)?\s*\*\s*Tie_min_neg\s*\)\s*\)",
        "Te_floor_exp", expression, flags=re.I,
    )
    expression = re.sub(
        r"\bdTe_floor_exp\b", "(Te_floor_exp/(0.5*Tie_min_neg))",
        expression, flags=re.I,
    )
    expression = re.sub(r"\bdTe_floor\b", "1", expression, flags=re.I)
    expression = re.sub(
        r"\bmin\s*\(\s*T0\s*,\s*T_min_neg\s*\)", "T0_floor",
        expression, flags=re.I,
    )
    expression = re.sub(
        r"\bexp\s*\(\s*\(\s*T0_floor\s*-\s*T_min_neg\s*\)\s*/\s*"
        r"\(\s*0\.5(?:d0)?\s*\*\s*T_min_neg\s*\)\s*\)",
        "T_floor_exp", expression, flags=re.I,
    )
    expression = re.sub(
        r"\bdT_floor_exp\b", "(T_floor_exp/(0.5*T_min_neg))",
        expression, flags=re.I,
    )
    expression = re.sub(r"\bdT_floor\b", "1", expression, flags=re.I)
    # Corrected neutral density, mirroring ``corr_neg_dens``.
    # ``rn0`` was renamed to ``rhon0`` at the top of this function, so the
    # corrected neutral density must be spelled the same way here; otherwise
    # the matcher compares ``rn0_corr`` with ``rhon0_corr`` while the report
    # displays both as ``rhon0_corr``.
    expression = re.sub(
        r"corr_neg_dens_n\(\s*rhon0?\s*\)", "rhon0_corr", expression,
        flags=re.I,
    )
    expression = re.sub(
        r"\bcorr_neg_dens_n\b", "rhon0_corr", expression, flags=re.I,
    )
    expression = re.sub(r"\bdrn0_corr_dn\b", "1", expression, flags=re.I)
    # Impurity negative-density correction, mirroring ``corr_neg_dens``.
    expression = re.sub(
        r"corr_neg_dens_imp\(\s*rhoimp0?\s*\)", "rhoimp0_corr",
        expression, flags=re.I,
    )
    expression = re.sub(
        r"\bcorr_neg_dens_imp\b", "rhoimp0_corr", expression, flags=re.I,
    )
    expression = re.sub(r"\bdrimp0_corr_dn\b", "1", expression, flags=re.I)
    expression = re.sub(
        r"\brimp0_corr\b", "rhoimp0_corr", expression, flags=re.I,
    )
    # ``sqrt`` is the only Fortran intrinsic appearing in the exported volume
    # terms.  Rewrite it as a power so the expression parser needs no function
    # table.
    expression = re.sub(
        r"\bsqrt\s*\(([^()]*)\)", r"((\1)**(1/2))", expression, flags=re.I,
    )
    # Fortran is case insensitive, and the routine spells the adiabatic index
    # both ``GAMMA`` and ``gamma`` — sometimes in a residual and its own
    # tangent.  Use one spelling.
    expression = re.sub(r"\bgamma\b", "GAMMA", expression, flags=re.I)
    # Fortran is case insensitive; the density and parallel-velocity blocks
    # spell the parallel-velocity trial function ``Vpar`` while the DSL prints
    # ``vpar``.  ``Vpar0`` is covered by its own alias below.
    expression = re.sub(
        r"\bVpar(_[a-z]+)?\b", r"vpar\1", expression, flags=re.I,
    )
    aliases = {
        "BB2": "((F0**2+ps0_x**2+ps0_y**2)/BigR**2)",
        "psi_grad2": "(ps0_x**2+ps0_y**2)",
        # Parallel-gradient work values of the density equation.  The density
        # ones are already covered by ``normalize_fortran_text``; the impurity
        # ones must be expanded after ``rimp0`` has been renamed.
        "Bgrad_rhoimp": "((F0*rhoimp0_p/BigR+rhoimp0_x*ps0_y-rhoimp0_y*ps0_x)/BigR)",
        "Bgrad_rhoimp_psi": "((rhoimp0_x*psi_y-rhoimp0_y*psi_x)/BigR)",
        "Bgrad_rhoimp_rhoimp": "((rhoimp_x*ps0_y-rhoimp_y*ps0_x)/BigR)",
        "Bgrad_rhoimp_rhoimp_n": "(F0*rhoimp_p/BigR**2)",
        # Parallel-velocity work values.
        # Single-temperature parallel-gradient work values.  ``Bgrad_T`` is
        # already defined by ``normalize_fortran_text``; only the tangents
        # are needed here.
        "Bgrad_T_T": "((T_x*ps0_y-T_y*ps0_x)/BigR)",
        "Bgrad_T_T_n": "(F0*T_p/BigR**2)",
        # Ion-temperature parallel-gradient work values.
        "Bgrad_Ti": "((F0*Ti0_p/BigR+Ti0_x*ps0_y-Ti0_y*ps0_x)/BigR)",
        "Bgrad_Ti_psi": "((Ti0_x*psi_y-Ti0_y*psi_x)/BigR)",
        "Bgrad_Ti_Ti": "((Ti_x*ps0_y-Ti_y*ps0_x)/BigR)",
        "Bgrad_Ti_Ti_n": "(F0*Ti_p/BigR**2)",
        # Electron-temperature parallel-gradient work values.
        "Bgrad_Te": "((F0*Te0_p/BigR+Te0_x*ps0_y-Te0_y*ps0_x)/BigR)",
        "Bgrad_Te_psi": "((Te0_x*psi_y-Te0_y*psi_x)/BigR)",
        "Bgrad_Te_Te": "((Te_x*ps0_y-Te_y*ps0_x)/BigR)",
        "Bgrad_Te_Te_n": "(F0*Te_p/BigR**2)",
        "Bgrad_vpar": "((F0*vpar0_p/BigR+vpar0_x*ps0_y-vpar0_y*ps0_x)/BigR)",
        "Bgrad_vpar_psi": "((vpar0_x*psi_y-vpar0_y*psi_x)/BigR)",
        "Bgrad_vpar_vpar": "((vpar_x*ps0_y-vpar_y*ps0_x)/BigR)",
        "Bgrad_vpar_vpar_n": "(F0*vpar_p/BigR**2)",
        # Prescribed rotation profile: a flux function whose flux derivative
        # is stored as ``dV_dpsi_source``.
        "Vt0_x": "(dV_dpsi_source*ps0_x)",
        "Vt0_y": "(dV_dpsi_source*ps0_y)",
        "Vt_x_psi": "(dV_dpsi_source*psi_x)",
        "Vt_y_psi": "(dV_dpsi_source*psi_y)",
        # Fortran is case insensitive: amat_n(var_vpar,var_Ti) spells the
        # background density ``R0``.
        "R0": "r0",
        "Btheta2": "((ps0_x**2+ps0_y**2)/BigR**2)",
        "Btheta2_psi": "(2*(psi_x*ps0_x+psi_y*ps0_y)/BigR**2)",
        "Vpar0": "vpar0",

        "Pe0": pe0,
        "Pe0_s": pe_derivatives["s"],
        "Pe0_t": pe_derivatives["t"],
        "Pe0_p": pe_derivatives["p"],
        "Pe0_x": pe_derivatives["x"],
        "Pe0_y": pe_derivatives["y"],
        "Pi0": pi0,
        "Pi0_s": pi_derivatives["s"],
        "Pi0_t": pi_derivatives["t"],
        "Pi0_y": pi_derivatives["y"],
        # Cartesian derivatives used by the diamagnetic viscosity block.
        # Keep these as explicit product rules so source and DSL expressions
        # are expanded to the same monomials by the report matcher.
        "Pi0_x": pi_derivatives["x"],
        "Pi0_xx": pi_derivatives["xx"],
        "Pi0_yy": pi_derivatives["yy"],
        "Pi0_xy": pi_derivatives["xy"],
        "P0": p0,
        "P0_s": p0_derivatives["s"],
        "P0_t": p0_derivatives["t"],
        "P0_p": p0_derivatives["p"],
    }
    for name, replacement in aliases.items():
        expression = re.sub(
            r"\b{}\b".format(name), replacement, expression, flags=re.I
        )
    if single_temperature:
        # The element routine sets Ti0 = Te0 = T0/2 in the one-temperature
        # branch (the evolved T0 is the total temperature), so a species
        # temperature reaching this point is half of the evolved one.
        expression = re.sub(r"\bTe0([_a-z]*)\b", r"(T0\1/2)", expression)
        expression = re.sub(r"\bTi0([_a-z]*)\b", r"(T0\1/2)", expression)
        # The one-temperature impurity closure stores the means of the
        # per-species coefficients:  alpha_imp = (alpha_i+alpha_e)/2,
        # alpha_imp_bis = (alpha_i+alpha_e_bis)/2 (alpha_i does not depend on
        # temperature, so it is its own "bis"), and alpha_imp_tri =
        # alpha_e_tri/4, the extra quarter coming from dTe0/dT0 = 1/2 applied
        # twice.  Express everything in the two-species basis so that
        # assignments shared between the branches and assignments written for
        # one branch only use the same symbols.
        expression = re.sub(
            r"\balpha_imp_tri\b", "(alpha_e_tri/4)", expression, flags=re.I,
        )
        expression = re.sub(
            r"\balpha_imp_bis\b", "((alpha_i+alpha_e_bis)/2)",
            expression, flags=re.I,
        )
        expression = re.sub(
            r"\balpha_imp\b", "((alpha_i+alpha_e)/2)", expression, flags=re.I,
        )
        # W_dia_Ti is assembled from Pi0_*_Ti, which the element routine
        # builds from the *trial* temperature basis function without the
        # dTi0/dT0 = 1/2 chain factor.  It is therefore twice the derivative
        # of W_dia with respect to the evolved one-temperature field.
        expression = re.sub(
            r"\bW_dia_Ti\b", "(2*W_dia_T)", expression, flags=re.I,
        )
    return expression
