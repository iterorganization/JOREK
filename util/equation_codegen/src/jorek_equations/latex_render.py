"""Render Python weak-form expressions (as ``ast`` nodes) in LaTeX.

The equations are written in Python (``model600.py``); this module turns their
source into display math for the documentation.  It does not import SymPy and
never evaluates the equations: it reads them as syntax, so a helper call such
as ``_B_dot_grad(rho)`` stays a named operator instead of being expanded.
"""

import ast
import copy
import json
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple


class RenderError(ValueError):
    """An expression that cannot be rendered."""

# ---------------------------------------------------------------------------
# LaTeX printing
# ---------------------------------------------------------------------------

ATOM, POW, MUL, ADD, LOW = 5, 4, 3, 2, 1

GREEK = {
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
    "kappa", "lambda", "mu", "nu", "xi", "pi", "rho", "sigma", "tau", "phi",
    "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi",
    "Sigma", "Phi", "Psi", "Omega",
}


def _name_part(part):
    match = re.fullmatch(r"([A-Za-z]+)(\d*)", part)
    if not match:
        return r"\mathrm{%s}" % part
    letters, digits = match.groups()
    if letters in GREEK:
        body = "\\" + letters
    elif len(letters) == 1:
        body = letters
    else:
        body = r"\mathrm{%s}" % letters
    return body + ("_{%s}" % digits if digits else "")


def default_symbol(name):
    """LaTeX for a name not in the notation table: ``D_par_local -> D_{par,local}``."""

    parts = [p for p in name.split("_") if p]
    if not parts:
        return r"\mathrm{%s}" % name.replace("_", r"\_")
    head = _name_part(parts[0])
    if len(parts) == 1:
        return head
    subscript = ",".join(_name_part(p) for p in parts[1:])
    if "_{" in head:
        return r"{%s}_{%s}" % (head, subscript)
    return r"%s_{%s}" % (head, subscript)


PLACEHOLDER = re.compile(r"\{(\w+)(!?)\}")

DERIVATIVES = {"dR": "R", "dZ": "Z", "dphi": r"\phi", "ds": "s", "dt": "t"}

BUILTIN_TEMPLATES = {
    "grad": (r"\nabla_{\mathrm{pol}} {0}", POW),
    "dot": (r"{0!}\cdot{1!}", MUL),
    "laplacian": (r"\nabla^2_{\mathrm{pol}} {0}", POW),
    "poiss_bracket": (r"[{0!},{1!}]", ATOM),
    # poiss_bracket_st(a, b) = J [a,b]; divided by xjac it is shown as [a,b]^{st}
    "poiss_bracket_st": (r"\mathcal{J}\,[{0!},{1!}]^{st}", MUL),
    "bracket": (r"\mathcal{J}\,[{0!},{1!}]", MUL),
    "__bracket_st": (r"[{0!},{1!}]^{st}", ATOM),
    "__bracket": (r"[{0!},{1!}]", ATOM),
    "freeze": (r"\overline{{0!}}", ATOM),
    "__dt": (r"\partial_t {0}", ATOM),
}


@dataclass
class Helper:
    name: str
    template: str
    params: List[str]
    keyword_only: List[str]
    defaults: Dict[str, ast.expr]
    display: Dict[str, str]


def _fill(template, values):
    """Substitute ``{key}``/``{key!}`` for the keys in ``values``; other braces are LaTeX."""

    def replace(match):
        key, raw = match.group(1), match.group(2)
        if key not in values:
            return match.group(0)
        latex, prec = values[key]
        if not raw and prec < ATOM:
            return r"\left(%s\right)" % latex
        return latex
    return PLACEHOLDER.sub(replace, template)


class LatexPrinter:
    """Render Python expressions and statements as LaTeX.

    ``values`` holds the local definitions of the function being rendered
    (name -> ``ast`` expression), used to see through ``psi_gradient[0]``;
    ``aliases`` maps a local ``functools.partial`` name to the helper it wraps.
    """

    def __init__(self, notation, helpers=None, local=None, terms=r"\Big(\textstyle\sum_k b_k\Big)",
                 values=None, aliases=None):
        self.notation = notation
        self.helpers = helpers or {}
        self.local = local or {}
        self.terms = terms
        self.values = values or {}
        self.aliases = aliases or {}

    def with_local(self, extra):
        return LatexPrinter(self.notation, self.helpers, dict(self.local, **extra), self.terms,
                            self.values, self.aliases)

    # -- names ---------------------------------------------------------------
    def symbol(self, name):
        if name in self.local:
            return self.local[name]
        if name in self.notation and "{0}" not in self.notation[name]:
            return self.notation[name]
        if name in BUILTIN_TEMPLATES:
            return _fill(BUILTIN_TEMPLATES[name][0], {"0": (r"\cdot", ATOM), "1": (r"\cdot", ATOM)})
        return default_symbol(name)

    # -- expressions -----------------------------------------------------------
    def expr(self, node) -> Tuple[str, int]:
        method = getattr(self, "_" + type(node).__name__, None)
        if method is None:
            raise RenderError("cannot render {} ({})".format(
                type(node).__name__, ast.unparse(node)))
        return method(node)

    def latex(self, node):
        return self.expr(node)[0]

    def wrap(self, node, minimum):
        latex, prec = self.expr(node)
        if prec < minimum:
            return r"\left(%s\right)" % latex
        return latex

    def _Constant(self, node):
        if isinstance(node.value, str):
            return r"\text{%s}" % node.value.replace("_", r"\_"), ATOM
        return repr(node.value), ATOM

    def _Name(self, node):
        if node.id == "terms":
            return self.terms, ATOM
        return self.symbol(node.id), ATOM

    def _Attribute(self, node):
        return self.symbol(node.attr), ATOM

    def _UnaryOp(self, node):
        sign = {ast.USub: "-", ast.UAdd: "+"}.get(type(node.op))
        if sign is None:
            raise RenderError("cannot render " + ast.unparse(node))
        return sign + self.wrap(node.operand, MUL), ADD

    def _BinOp(self, node):
        op = type(node.op)
        if op is ast.Add:
            return "%s + %s" % (self.wrap(node.left, ADD), self.wrap(node.right, MUL if _is_signed(node.right) else ADD)), ADD
        if op is ast.Sub:
            return "%s - %s" % (self.wrap(node.left, ADD), self.wrap(node.right, MUL)), ADD
        if op in (ast.Mult, ast.Div):
            sign, unsigned = _pull_sign(node)
            if sign == "-":
                return "-" + self.wrap(unsigned, MUL), ADD
            node = unsigned
        if op in (ast.Mult, ast.Div):
            numerator, denominator = _cancel_jacobian(*_num_den(node))
            if not denominator:
                return self.product(numerator), MUL
            bottom = self.latex(denominator[0]) if len(denominator) == 1 else self.product(denominator)
            number = 1
            for factor in numerator:
                if _is_number(factor):
                    number *= factor.value
            rest = [f for f in numerator if not _is_number(f)]
            if len(rest) == 1 and self.expr(rest[0])[1] <= ADD:
                return r"\frac{%s}{%s}\left(%s\right)" % (repr(number), bottom, self.latex(rest[0])), MUL
            if len(rest) <= 1:
                top = self.latex(rest[0]) if rest and number == 1 else self.product(numerator)
                return r"\frac{%s}{%s}" % (top, bottom), MUL
            # keep the factors at full size: (number / denominator) * factors
            return r"\frac{%s}{%s}\,%s" % (repr(number), bottom, self.product(rest)), MUL
        if op is ast.Pow:
            return "%s^{%s}" % (self.wrap(node.left, ATOM), self.latex(node.right)), POW
        raise RenderError("cannot render " + ast.unparse(node))

    def product(self, factors):
        """Juxtapose factors, numbers first; bracket dot products that are followed by more."""

        number = 1
        for factor in factors:
            if _is_number(factor):
                number *= factor.value
        rest = [f for f in factors if not _is_number(f)]
        rest = [f for f in rest if _is_rational(f)] + [f for f in rest if not _is_rational(f)]
        parts = [repr(number)] if number != 1 or not rest else []
        for position, factor in enumerate(rest):
            latex = self.wrap(factor, MUL)
            if (re.search(r"\\cdot(?![a-z])", latex) and len(rest) > 1
                    and not latex.startswith(r"\left(")):
                latex = r"\left(%s\right)" % latex
            parts.append(latex)
        return r"\,".join(parts)

    def _IfExp(self, node):
        return (
            r"\begin{cases} %s & \text{if } %s \\ %s & \text{otherwise} \end{cases}"
            % (self.latex(node.body), _condition(node.test), self.latex(node.orelse))
        ), ATOM

    def _Lambda(self, node):
        params = [a.arg for a in node.args.args]
        inner = self.with_local({p: p for p in params})
        return r"(%s) \mapsto %s" % (", ".join(params), inner.latex(node.body)), LOW

    def _Tuple(self, node):
        return r"\left(%s\right)" % ", ".join(self.latex(e) for e in node.elts), ATOM

    def _Subscript(self, node):
        value, index = node.value, node.slice
        if isinstance(value, ast.Name) and isinstance(self.values.get(value.id), ast.Call):
            value = self.values[value.id]
        if (isinstance(value, ast.Call) and _call_name(value.func) == "grad" and len(value.args) == 1
                and isinstance(index, ast.Constant) and index.value in (0, 1)):
            return r"\partial_{%s} %s" % ("RZ"[index.value], self.wrap(value.args[0], ATOM)), POW
        return "{%s}_{%s}" % (self.wrap(node.value, ATOM), self.latex(node.slice)), ATOM

    def _GeneratorExp(self, node):
        """``(f(a, b) for a, b in zip(A, B))``: the vector ``f(A, B)``."""

        (loop,) = node.generators
        if (loop.ifs or not isinstance(loop.iter, ast.Call) or _call_name(loop.iter.func) != "zip"
                or not isinstance(loop.target, ast.Tuple)):
            raise RenderError("cannot render " + ast.unparse(node))
        names = [target.id for target in loop.target.elts]
        inner = self.with_local({n: self.latex(a) for n, a in zip(names, loop.iter.args)})
        return inner.expr(node.elt)

    def _Call(self, node):
        name = _call_name(node.func)
        name = self.aliases.get(name, name)
        args = node.args
        if name == "tuple" and len(args) == 1 and isinstance(args[0], ast.GeneratorExp):
            return self.expr(args[0])
        if name in DERIVATIVES and len(args) == 1:
            chain = [DERIVATIVES[name]]
            inner = args[0]
            while (isinstance(inner, ast.Call) and _call_name(inner.func) in DERIVATIVES
                   and len(inner.args) == 1):
                chain.append(DERIVATIVES[_call_name(inner.func)])
                inner = inner.args[0]
            ops = "".join(r"\partial_{%s}" % c for c in chain)
            return "%s %s" % (ops, self.wrap(inner, ATOM)), POW
        if name == "Rational" and len(args) == 2:
            return r"\frac{%s}{%s}" % (self.latex(args[0]), self.latex(args[1])), ATOM
        if name == "sqrt" and len(args) == 1:
            return r"\sqrt{%s}" % self.latex(args[0]), ATOM
        if name in ("coefficient", "test_function") and args and isinstance(args[0], ast.Constant):
            return self.symbol(args[0].value), ATOM
        if name in BUILTIN_TEMPLATES:
            template, prec = BUILTIN_TEMPLATES[name]
            return _fill(template, {str(i): self.expr(a) for i, a in enumerate(args)}), prec
        if name in self.helpers:
            helper = self.helpers[name]
            names = set(helper.params + helper.keyword_only)
            shown = [k for k, _ in PLACEHOLDER.findall(helper.template) if k in names]
            atomic = not shown or helper.template.rstrip().endswith(")")
            return self._helper_call(helper, node), ATOM if atomic else POW
        values = [self.latex(a) for a in args]
        if name in self.notation and "{0}" in self.notation[name]:
            return _fill(self.notation[name], {str(i): self.expr(a) for i, a in enumerate(args)}), POW
        head = self.symbol(name)
        if not values:
            return head, ATOM
        return r"%s(%s)" % (head, ", ".join(values)), ATOM

    def _helper_call(self, helper, node):
        values = {}
        for param, arg in zip(helper.params, node.args):
            values[param] = arg
        for keyword in node.keywords:
            values[keyword.arg] = keyword.value
        names = set(helper.params + helper.keyword_only)
        used = {key for key, _ in PLACEHOLDER.findall(helper.template) if key in names}
        for param in values:
            if param not in used and param not in helper.keyword_only:
                raise RenderError(
                    "{}(...) is called with {!r}, which its template {!r} does not show; "
                    "add {{{}}} to the template".format(helper.name, param, helper.template, param))
        filled = {}
        for key in used:
            if key in values:
                filled[key] = self.expr(values[key])
            elif key in helper.defaults:
                filled[key] = self.expr(helper.defaults[key])
            else:
                raise RenderError("{}(...) is called without {!r}".format(helper.name, key))
        return _fill(helper.template, filled)

    # -- statements ------------------------------------------------------------
    def rows(self, statements, indent=""):
        rows = []
        for statement in statements:
            if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Constant):
                continue
            if isinstance(statement, ast.Raise) or (
                    isinstance(statement, ast.If) and not statement.orelse
                    and all(isinstance(s, ast.Raise) for s in statement.body)):
                continue  # guards are not part of the equation
            if (isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Call)
                    and _call_name(statement.value.func) == "partial"):
                continue  # an alias of a helper, resolved where it is called
            if isinstance(statement, ast.Assign):
                targets = " = ".join(self.latex(t) for t in statement.targets)
                value = self.latex(statement.value)
                if value != targets:
                    rows.append(r"%s%s &= %s" % (indent, targets, value))
            elif isinstance(statement, ast.AugAssign) and isinstance(statement.op, (ast.Add, ast.Sub)):
                op = "+" if isinstance(statement.op, ast.Add) else "-"
                rows.append(r"%s%s &\mathrel{%s}= %s" % (
                    indent, self.latex(statement.target), op, self.latex(statement.value)))
            elif isinstance(statement, ast.If):
                rows.append(r"&%s\text{if } %s\text{:}" % (indent, _condition(statement.test)))
                rows.extend(self.rows(statement.body, indent + r"\quad "))
                if statement.orelse:
                    rows.append(r"&%s\text{otherwise:}" % indent)
                    rows.extend(self.rows(statement.orelse, indent + r"\quad "))
            elif isinstance(statement, ast.Return):
                rows.append(r"%s &= %s" % (indent, self.latex(statement.value)))
            else:
                raise RenderError("cannot render statement " + ast.unparse(statement))
        return rows


def _call_name(func):
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return ""


def _is_signed(node):
    return isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd))


def _is_number(node):
    return isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) \
        and not isinstance(node.value, bool)


def _is_rational(node):
    return isinstance(node, ast.Call) and _call_name(node.func) == "Rational"


def _pull_sign(node):
    """Move a unary sign on the leftmost factor of a product to the front."""

    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        return ("-" if isinstance(node.op, ast.USub) else "+"), node.operand
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Mult, ast.Div)):
        sign, left = _pull_sign(node.left)
        if left is not node.left:
            return sign, ast.BinOp(left=left, op=node.op, right=node.right)
    return "+", node


def _num_den(node):
    """Numerator and denominator factors of a product/quotient chain."""

    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        left, right = _num_den(node.left), _num_den(node.right)
        return left[0] + right[0], left[1] + right[1]
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left, right = _num_den(node.left), _num_den(node.right)
        return left[0] + right[1], left[1] + right[0]
    return [node], []


BRACKETS_ST = {"poiss_bracket_st": "__bracket_st", "bracket": "__bracket"}


def _cancel_jacobian(numerator, denominator):
    """Pair each element-coordinate bracket with a ``xjac`` of the denominator."""

    numerator, denominator = list(numerator), list(denominator)
    if any(_call_name(getattr(f, "func", None)) in BRACKETS_ST for f in numerator if isinstance(f, ast.Call)):
        # the volume weight dV = R*xjac in a denominator also supplies a xjac
        expanded = []
        for factor in denominator:
            if isinstance(factor, ast.Name) and factor.id == "dV":
                expanded += [ast.Name(id="R", ctx=ast.Load()), ast.Name(id="xjac", ctx=ast.Load())]
            else:
                expanded.append(factor)
        denominator = expanded
    for index, factor in enumerate(numerator):
        name = _call_name(factor.func) if isinstance(factor, ast.Call) else ""
        if name not in BRACKETS_ST:
            continue
        jacobian = next((i for i, d in enumerate(denominator)
                         if isinstance(d, ast.Name) and d.id == "xjac"), None)
        if jacobian is None:
            continue
        del denominator[jacobian]
        numerator[index] = ast.Call(func=ast.Name(id=BRACKETS_ST[name], ctx=ast.Load()),
                                    args=factor.args, keywords=[])
    return numerator, denominator


def _rebuild(numerator, denominator):
    def chain(factors):
        node = factors[0]
        for factor in factors[1:]:
            node = ast.BinOp(left=node, op=ast.Mult(), right=factor)
        return node
    top = chain(numerator) if numerator else ast.Constant(value=1)
    return ast.BinOp(left=top, op=ast.Div(), right=chain(denominator)) if denominator else top


def _without_volume_element(node, name="dV"):
    """``node / dV`` for an expression with a ``dV`` factor in its numerator."""

    numerator, denominator = _num_den(node)
    for index, factor in enumerate(numerator):
        if isinstance(factor, ast.Name) and factor.id == name:
            return _rebuild(numerator[:index] + numerator[index + 1:], denominator)
    raise RenderError("expected a factor {} in {}".format(name, ast.unparse(node)))


def _factors(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        return _factors(node.left) + _factors(node.right)
    return [node]


def _condition(test):
    # GitHub's Markdown-to-KaTeX pipeline unescapes "\_" back to "_" before the
    # math is parsed, which then errors inside \texttt{...} (text mode, where a
    # bare "_" is invalid).  Render the flag name with spaces instead of a
    # backslash escape, so no "\_" sequence reaches the page at all.
    return r"\texttt{%s}" % ast.unparse(test).replace("_", " ")


def _signed_terms(node):
    """Split a sum into ``(sign, term)`` pairs, keeping the written order."""

    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        left = _signed_terms(node.left)
        right = _signed_terms(node.right)
        if isinstance(node.op, ast.Sub):
            right = [("-" if s == "+" else "+", t) for s, t in right]
        return left + right
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        sign = "-" if isinstance(node.op, ast.USub) else "+"
        return [(sign if s == "+" else ("+" if sign == "-" else "-"), t)
                for s, t in _signed_terms(node.operand)] if _is_sum(node.operand) else [(sign, node.operand)]
    sign, unsigned = _pull_sign(node)
    if unsigned is not node:
        return [(sign, unsigned)]
    return [("+", node)]


def _is_sum(node):
    return isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub))


def _aligned(rows):
    if len(rows) == 1 and "&" not in rows[0].replace(r"\&", ""):
        return rows[0]
    return "\\begin{aligned}\n" + " \\\\\n".join(rows) + "\n\\end{aligned}"


# ---------------------------------------------------------------------------
# Helpers and blocks
# ---------------------------------------------------------------------------


LONG_SUM = 4
LONG_LINE = 160  # characters of LaTeX that still read well on one line


def _definition_rows(printer, lhs, value):
    """``lhs = value``, one term per row when the value is a long sum."""

    single = r"%s &= %s" % (lhs, printer.latex(value))
    if len(single) <= LONG_LINE:
        return [single]
    terms = _signed_terms(value)
    if len(terms) >= LONG_SUM:
        rows = [r"%s &= %s" % (lhs, _signed(printer, terms[0], first=True))]
        rows += [r"&\quad %s" % _signed(printer, term) for term in terms[1:]]
        return rows
    sign, unsigned = _pull_sign(value)
    numerator, denominator = _num_den(unsigned)
    if not denominator and len(numerator) > 1:
        inner = _signed_terms(numerator[-1])
        if len(inner) >= 3:
            factor = ("-" if sign == "-" else "") + printer.product(numerator[:-1])
            rows = [r"%s &= %s\Big(%s" % (lhs, factor, _signed(printer, inner[0], first=True))]
            rows += [r"&\qquad %s" % _signed(printer, term) for term in inner[1:]]
            rows[-1] += r"\Big)"
            return rows
    return [r"%s &= %s" % (lhs, printer.latex(value))]


def _signed(printer, signed_term, first=False):
    sign, term = signed_term
    latex = printer.wrap(term, MUL)
    if first:
        return ("-" if sign == "-" else "") + latex
    return "%s %s" % (sign, latex)


def _is_variant(statement):
    """An ``if flag: return ...`` spelling variant of a helper's result."""

    return (isinstance(statement, ast.If) and not statement.orelse
            and all(isinstance(s, ast.Return) for s in statement.body))


SENTINEL = "\x00terms\x00"
FIELD_NAMES = {"psi", "u", "j", "omega", "rho", "T", "Ti", "Te", "vpar", "rhon", "rhoimp"}


class _TimeDerivative:
    """Write out the time derivative of a mass term.

    ``freeze(x)`` marks a factor that the time derivative does not act on, so
    a mass term with frozen factors is expanded by the product rule over its
    remaining factors.  A term without frozen factors is kept as
    ``partial_t(term)``.
    """

    def __init__(self, values, helpers):
        self.values = values
        self.helpers = helpers

    def depends(self, node, seen=()):
        if isinstance(node, ast.Call) and _call_name(node.func) == "freeze":
            return False
        if isinstance(node, ast.Name):
            if node.id in FIELD_NAMES:
                return True
            if node.id in self.values and node.id not in seen:
                return self.depends(self.values[node.id], seen + (node.id,))
            return False
        if isinstance(node, ast.Call) and _call_name(node.func) in self.helpers:
            helper = self.helpers[_call_name(node.func)]
            given = set(helper.params[:len(node.args)]) | {k.arg for k in node.keywords}
            if any(self.depends(d, seen) for p, d in helper.defaults.items()
                   if p not in given and p not in helper.keyword_only):
                return True
        return any(self.depends(child, seen) for child in ast.iter_child_nodes(node))

    @staticmethod
    def has_freeze(node):
        return any(isinstance(n, ast.Call) and _call_name(n.func) == "freeze" for n in ast.walk(node))

    @staticmethod
    def thaw(node):
        class Thaw(ast.NodeTransformer):
            def visit_Call(self, call):
                self.generic_visit(call)
                if _call_name(call.func) == "freeze" and len(call.args) == 1:
                    return call.args[0]
                return call
        import copy
        return Thaw().visit(copy.deepcopy(node))

    LINEAR = {"grad", "dR", "dZ", "dphi", "laplacian"}

    def mark(self, node):
        """``partial_t node``, moved inside linear operators of a single time-dependent argument."""

        if isinstance(node, ast.Name) and node.id in self.values and node.id not in FIELD_NAMES:
            value = self.values[node.id]
            if isinstance(value, ast.Call) and _call_name(value.func) in self.LINEAR | {"dot"}:
                return self.mark(value)
        if isinstance(node, ast.Call) and not node.keywords:
            name = _call_name(node.func)
            timed = [i for i, a in enumerate(node.args) if self.depends(a)]
            if (name in self.LINEAR or name == "dot") and len(timed) == 1:
                args = list(node.args)
                args[timed[0]] = self.mark(args[timed[0]])
                return ast.Call(func=node.func, args=args, keywords=[])
        return ast.Call(func=ast.Name(id="__dt", ctx=ast.Load()), args=[node], keywords=[])

    def derive(self, node):
        if not self.depends(node):
            return None
        if not self.has_freeze(node):
            return self.mark(node)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
            left, right = self.derive(node.left), self.derive(node.right)
            if right is None:
                return left
            if left is None:
                return right if isinstance(node.op, ast.Add) else ast.UnaryOp(op=ast.USub(), operand=right)
            return ast.BinOp(left=left, op=node.op, right=right)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            inner = self.derive(node.operand)
            return None if inner is None else ast.UnaryOp(op=node.op, operand=inner)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Mult, ast.Div)):
            numerator, denominator = _num_den(node)
            if any(self.depends(d) for d in denominator):
                return self.mark(self.thaw(node))
            result = None
            for index, factor in enumerate(numerator):
                derived = self.derive(factor)
                if derived is None:
                    continue
                others = [self.thaw(f) for f in numerator[:index] + numerator[index + 1:]]
                term = _rebuild(others + [derived], denominator)
                result = term if result is None else ast.BinOp(left=result, op=ast.Add(), right=term)
            return result
        return self.mark(self.thaw(node))



POLICY_TEXT = {
    "active": "differentiated",
    "piecewise_active": "differentiated (current branch)",
    "frozen": "not differentiated",
}


def functions_table(symbols_source):
    """A Markdown table of the ``external_function`` declarations of a module."""

    rows = ["| DSL | Arguments | Fortran value | Derivatives supplied | Jacobian |",
            "|---|---|---|---|---|"]
    for statement in ast.parse(symbols_source).body:
        if not (isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Call)
                and _call_name(statement.value.func) == "external_function"):
            continue
        call = statement.value
        keywords = {k.arg: ast.literal_eval(k.value) for k in call.keywords}
        name = statement.targets[0].id
        fortran = keywords.get("fortran_name") or ast.literal_eval(call.args[0])
        arguments = ", ".join("`%s`" % a for a in keywords.get("arguments", ())) or "none"
        derivatives = ", ".join("`%s` for `%s`" % (d, a) for a, d in keywords.get("derivatives", {}).items())
        policy = POLICY_TEXT[keywords.get("policy", "active")]
        rows.append("| `%s` | %s | `%s` | %s | %s |" % (name, arguments, fortran, derivatives or "none", policy))
    return rows


FINDINGS_URL = ("https://github.com/iterorganization/JOREK/blob/develop/"
                "util/equation_codegen/JOREK_FINDINGS.md")


def findings_table(reference_source):
    """A Markdown table of the findings recorded in the checker's reference."""

    import json

    blocks = json.loads(reference_source)["blocks"]
    by_finding = {}
    for name, block in blocks.items():
        lhs = name.split(" | ")[1]
        for finding in block["findings"]:
            kinds = by_finding.setdefault(finding, ({}, {}))
            target = kinds[0] if lhs.startswith("rhs") else kinds[1]
            target[lhs] = target.get(lhs, 0) + len(block["source_only"]) + len(block["generated_only"])
    rows = ["| Finding | Residual blocks | Jacobian blocks |", "|---:|---|---|"]
    for finding in sorted(by_finding):
        residual, jacobian = by_finding[finding]
        rows.append("| [%d](%s) | %s | %s |" % (
            finding, FINDINGS_URL,
            ", ".join("`%s`" % b for b in sorted(residual)) or "none",
            ", ".join("`%s`" % b for b in sorted(jacobian)) or "none"))
    return rows

