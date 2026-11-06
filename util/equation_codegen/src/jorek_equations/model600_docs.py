"""Generate the model-600 weak-form documentation from ``model600.py``.

``model600.py`` is the only place the equations are written.  This module
reads its source (with ``ast``, without importing or evaluating it) and writes
the generated regions of ``docs/physics/base_fluid_models/RMHD/weak_form.md``:
the full equations, the notation tables, the shared operators, the
per-equation term groups, and two tables built from the checker's data.  The
prose around them is hand-written and kept as it is.

A generated region is delimited by::

    <!-- BEGIN GENERATED: name -->
    ...
    <!-- END GENERATED: name -->
"""

import ast
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import model600_notation as notation
from .latex_render import (
    ATOM,
    MUL,
    SENTINEL,
    Helper,
    LatexPrinter,
    RenderError,
    _TimeDerivative,
    _aligned,
    _call_name,
    _definition_rows,
    _fill,
    _is_variant,
    _num_den,
    _rebuild,
    _signed_terms,
    functions_table,
    findings_table,
)

SOURCE_NAME = "util/equation_codegen/src/jorek_equations/model600.py"
BEGIN = re.compile(r"^<!-- BEGIN GENERATED: ([\w-]+) -->\s*$")
END = "<!-- END GENERATED: {} -->"
MATH = "$$"


# ---------------------------------------------------------------------------
# Reading model600.py
# ---------------------------------------------------------------------------

@dataclass
class TermGroup:
    title: str
    terms: List[Tuple[str, ast.expr]]


@dataclass
class EquationSource:
    name: str
    title: str
    anchor: str
    prose: str
    definitions: List[ast.stmt]
    values: Dict[str, ast.expr]
    aliases: Dict[str, str]
    mass: Optional[ast.expr]
    rhs_name: str
    rhs: Optional[ast.expr]
    groups: List[TermGroup] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def _title(text):
    text = text.strip().rstrip(".").strip()
    return text[:1].upper() + text[1:]


def _comment_block(lines, lineno):
    """The comment lines directly above line ``lineno`` (1-based), top first."""

    comments, index = [], lineno - 2
    while index >= 0 and lines[index].strip().startswith("#"):
        comments.insert(0, lines[index].strip().lstrip("#").strip())
        index -= 1
    return comments


def _split_groups(inner_lines):
    """Split the lines inside ``B = ( ... )`` into ``(title, code)`` groups at comments."""

    groups, title, code, depth = [], None, [], 0
    for line in inner_lines:
        stripped = line.strip()
        if depth == 0 and stripped.startswith("#"):
            comment = stripped.lstrip("#").strip()
            if code:
                groups.append((title, code))
                title, code = comment, []
            else:
                title = comment if title is None else title + " " + comment
            continue
        if depth == 0 and not stripped:
            if code:
                groups.append((title, code))
                title, code = None, []
            continue
        code.append(line)
        cleaned = re.sub(r"#.*", "", line)
        depth += cleaned.count("(") - cleaned.count(")")
    if code:
        groups.append((title, code))
    return [(_title(t) if t else "Other terms", "\n".join(c)) for t, c in groups]


def _local_names(statements):
    values, aliases = {}, {}
    for statement in ast.walk(ast.Module(body=list(statements), type_ignores=[])):
        if isinstance(statement, ast.Assign) and len(statement.targets) == 1 \
                and isinstance(statement.targets[0], ast.Name):
            name, value = statement.targets[0].id, statement.value
            if isinstance(value, ast.Call) and _call_name(value.func) == "partial" and value.args:
                aliases[name] = _call_name(value.args[0])
            else:
                values.setdefault(name, value)
    return values, aliases


def read_equation(function, lines, title, anchor, prose):
    body = function.body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        body = body[1:]
    values, aliases = _local_names(body)
    equation = EquationSource(
        name=function.name, title=title, anchor=anchor, prose=prose,
        definitions=[], values=values, aliases=aliases, mass=None, rhs_name="", rhs=None,
    )
    for statement in body:
        target = (statement.targets[0].id if isinstance(statement, ast.Assign)
                  and isinstance(statement.targets[0], ast.Name) else None)
        if isinstance(statement, ast.Return):
            equation.notes = _comment_block(lines, _first_code_line_above(lines, statement.lineno))
            continue
        if target == "A":
            equation.mass = statement.value
        elif target in ("B", "C"):
            equation.rhs_name = target
            equation.rhs = statement.value
            text = lines[statement.lineno - 1:statement.end_lineno]
            if target == "B":
                if text[0].strip() != "B = (" or not text[-1].strip().startswith(")"):
                    raise RenderError("{}: write the right-hand side as 'B = (' ... ') * dV'".format(
                        function.name))
                for title_, code in _split_groups(text[1:-1]):
                    expression = ast.parse("(\n" + _dedent(code) + "\n)", mode="eval").body
                    equation.groups.append(TermGroup(title_, _signed_terms(expression)))
        else:
            equation.definitions.append(statement)
    if equation.rhs is None:
        raise RenderError("{} defines neither B nor C".format(function.name))
    return equation


def _first_code_line_above(lines, lineno):
    """Skip the blank lines above ``lineno`` so a comment block separated by one is found."""

    while lineno > 1 and not lines[lineno - 2].strip():
        lineno -= 1
    return lineno


def _dedent(code):
    import textwrap
    return textwrap.dedent(code).strip("\n")


def read_helpers(tree):
    helpers, functions = {}, {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in notation.HELPERS:
            info = notation.HELPERS[node.name]
            arguments = node.args
            params = [a.arg for a in arguments.args]
            defaults = dict(zip(params[len(params) - len(arguments.defaults):], arguments.defaults))
            for arg, default in zip(arguments.kwonlyargs, arguments.kw_defaults):
                if default is not None:
                    defaults[arg.arg] = default
            helpers[node.name] = Helper(
                name=node.name, template=info["template"], params=params,
                keyword_only=[a.arg for a in arguments.kwonlyargs],
                defaults=defaults, display=dict(info["display"]),
            )
            functions[node.name] = node
    missing = [name for name in notation.HELPERS if name not in helpers]
    if missing:
        raise RenderError("helpers documented but not in model600.py: " + ", ".join(missing))
    undocumented = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name.startswith("_") and n.name not in notation.HELPERS]
    if undocumented:
        raise RenderError("add these helpers to model600_notation.HELPERS: " + ", ".join(undocumented))
    return helpers, functions


def read_source(source):
    tree = ast.parse(source)
    lines = source.split("\n")
    helpers, helper_functions = read_helpers(tree)
    functions = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    equations = []
    for name, title, anchor, prose in notation.EQUATIONS:
        if name not in functions:
            raise RenderError("equation {} is documented but not in model600.py".format(name))
        equations.append(read_equation(functions[name], lines, title, anchor, prose))
    public = [n for n in functions if not n.startswith("_")]
    missing = [n for n in public if n not in {e.name for e in equations}]
    if missing:
        raise RenderError("add these equations to model600_notation.EQUATIONS: " + ", ".join(missing))
    return helpers, helper_functions, equations


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def symbol_map():
    table = {}
    table.update(notation.FIELDS)
    table.update(notation.SYMBOLS)
    table.update(notation.FUNCTIONS)
    return table


def _without_weight(node):
    """``node / (R*xjac)``: drop the volume weight ``dV = R*xjac`` from a term."""

    numerator, denominator = _num_den(node)
    names = [f.id if isinstance(f, ast.Name) else None for f in numerator]
    if "dV" in names:
        del numerator[names.index("dV")]
        return _rebuild(numerator, denominator)
    if "xjac" not in names:
        raise RenderError("no volume weight in " + ast.unparse(node))
    del numerator[names.index("xjac")]
    names = [f.id if isinstance(f, ast.Name) else None for f in numerator]
    if "R" in names:
        del numerator[names.index("R")]
    else:
        below = [f.id if isinstance(f, ast.Name) else None for f in denominator]
        if "R" in below:
            denominator[below.index("R")] = ast.BinOp(
                left=ast.Name(id="R", ctx=ast.Load()), op=ast.Pow(), right=ast.Constant(value=2))
        else:
            denominator.append(ast.Name(id="R", ctx=ast.Load()))
    return _rebuild(numerator, denominator)


def _math(latex):
    return [MATH] + latex.split("\n") + [MATH]


def render_summary(equation, symbols, helpers):
    printer = LatexPrinter(symbols, helpers, values=equation.values, aliases=equation.aliases)
    if equation.rhs_name == "B":
        rows = ["%s %s" % ("" if i == 0 and sign == "+" else sign, printer.wrap(term, MUL))
                for i, (sign, term) in enumerate(t for g in equation.groups for t in g.terms)]
        mass = _without_weight(equation.mass)
        if _TimeDerivative.has_freeze(mass):
            lhs = printer.latex(_TimeDerivative(equation.values, helpers).derive(mass))
        else:
            lhs = r"\frac{\partial}{\partial t}\Big(%s\Big)" % printer.latex(mass)
    else:
        lhs, rows = "0", [printer.latex(_without_weight(equation.rhs))]
    rows[0] = r"%s \;=\; &%s" % (lhs, rows[0])
    rows[1:] = ["&" + row for row in rows[1:]]
    return "\\begin{aligned}\n" + " \\\\\n".join(rows) + "\n\\end{aligned}"


def render_helper(function, helper, symbols, helpers):
    printer = LatexPrinter(symbols, helpers)
    local = {}
    for param in helper.params + helper.keyword_only:
        if param in helper.display:
            local[param] = helper.display[param]
        elif param in helper.defaults and param not in helper.keyword_only:
            local[param] = printer.latex(helper.defaults[param])
    inner = printer.with_local(local)
    lhs = _fill(helper.template, {p: (local.get(p, inner.symbol(p)), ATOM)
                                  for p in helper.params + helper.keyword_only})
    body = function.body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        body = body[1:]
    body = [s for s in body if not _is_variant(s)]
    if not isinstance(body[-1], ast.Return):
        raise RenderError("helper {} must end with a return".format(function.name))
    rows = inner.rows(body[:-1])
    rows.extend(_definition_rows(inner, lhs, body[-1].value))
    return _aligned(rows)


def region_equations(equations, symbols, helpers):
    out = []
    for equation in equations:
        out += ['<a id="eq-%s"></a>' % equation.anchor, "", "### " + equation.title, ""]
        out += _math(render_summary(equation, symbols, helpers))
        out += ["", "[Definitions and term groups](#details-%s)" % equation.anchor, ""]
    return out


def region_operators(helper_functions, symbols, helpers):
    out, group = [], None
    for name, info in notation.HELPERS.items():
        if info["group"] != group:
            group = info["group"]
            out += ["### " + group, ""]
        out += ["#### " + info["title"], ""]
        if info["description"]:
            out += [info["description"], ""]
        out += _math(render_helper(helper_functions[name], helpers[name], symbols, helpers))
        out += ["", "Source: `%s`." % name, ""]
    return out


def region_details(equations, symbols, helpers):
    out = []
    for equation in equations:
        printer = LatexPrinter(symbols, helpers, values=equation.values, aliases=equation.aliases)
        out += ['<a id="details-%s"></a>' % equation.anchor, "", "### " + equation.title, ""]
        if equation.prose:
            out += [equation.prose, ""]
        out += ["Source: `%s` in `%s`." % (equation.name, SOURCE_NAME), ""]
        definitions = printer.rows(equation.definitions)
        if definitions:
            out += ["#### Definitions", ""] + _math(_aligned(definitions)) + [""]
        if equation.mass is not None:
            out += ["#### Time derivative", ""]
            out += _math(r"A = %s" % printer.latex(equation.mass)) + [""]
        if equation.rhs_name == "C":
            out += ["#### Constraint", ""] + _math(r"C = %s" % printer.latex(equation.rhs)) + [""]
            continue
        out += ["#### Right-hand side", ""]
        out += _math(r"B = %s" % LatexPrinter(symbols, helpers).latex(
            _replace_sum(equation.rhs))) + [""]
        for group in equation.groups:
            out += ["##### " + group.title, ""]
            rows = ["&%s %s" % (sign, printer.wrap(term, MUL)) for sign, term in group.terms]
            out += _math("\\begin{aligned}\n" + " \\\\\n".join(rows) + "\n\\end{aligned}") + [""]
        if equation.notes:
            out += ["**Note:** " + " ".join(equation.notes), ""]
    return out


def _replace_sum(node):
    """``B = (sum) * dV`` shown as ``(sum_k b_k) dV``."""

    return ast.BinOp(left=ast.Name(id="terms", ctx=ast.Load()), op=node.op, right=node.right)


def region_notation():
    out = ["### Work values", "", "| Python | Symbol |", "|---|---|"]
    out += ["| `%s` | $%s$ |" % (n, s) for n, s in notation.SYMBOLS.items()]
    out += ["", "### State-dependent functions", "",
            r"A function of the state is shown with its arguments, e.g. $S_{ion}(T)$;",
            r"$x^{\mathrm{c}}$ is a density corrected for negative values.", "",
            "| Python | Symbol |", "|---|---|"]
    out += ["| `%s` | $%s$ |" % (n, s) for n, s in notation.FUNCTIONS.items()]
    return out


def render_regions(source, reference):
    helpers, helper_functions, equations = read_source(source)
    symbols = symbol_map()
    return {
        "equations": region_equations(equations, symbols, helpers),
        "notation": region_notation(),
        "linearization": functions_table(source),
        "findings": findings_table(reference),
        "operators": region_operators(helper_functions, symbols, helpers),
        "details": region_details(equations, symbols, helpers),
    }


def render_page(text, source, reference):
    """Return ``text`` with every generated region rewritten."""

    regions = render_regions(source, reference)
    lines = text.split("\n")
    output, index, seen = [], 0, set()
    while index < len(lines):
        match = BEGIN.match(lines[index])
        if not match:
            output.append(lines[index])
            index += 1
            continue
        name = match.group(1)
        if name not in regions:
            raise RenderError("unknown generated region {!r}".format(name))
        end = index + 1
        while end < len(lines) and lines[end].strip() != END.format(name):
            end += 1
        if end == len(lines):
            raise RenderError("generated region {!r} is not closed".format(name))
        output += [lines[index], ""] + regions[name] + ["", lines[end]]
        seen.add(name)
        index = end + 1
    missing = set(regions) - seen
    if missing:
        raise RenderError("the page has no region for: " + ", ".join(sorted(missing)))
    return "\n".join(output)
