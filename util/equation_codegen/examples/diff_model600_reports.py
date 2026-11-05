#!/usr/bin/env python3
"""Compare the two model-600 term reports block by block.

Meld shows every line that differs, including the many lines that differ only
because source and generated monomials are listed in a different order.  This
script instead compares each assignment block as a multiset of monomials and,
where the multisets differ, prints the algebraic residual ``source-generated``.
A block whose residual is zero contains the same physics in a different
arrangement; a block with a non-zero residual is a real disagreement.
"""

import collections
from collections import Counter, OrderedDict
from pathlib import Path
import argparse
import re
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def read_blocks(path):
    """Return {(section, assignment, repeat): [monomial, ...]}."""

    blocks = OrderedDict()
    key = None
    section = ""
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or line.startswith(">"):
            # Any heading closes the block being read.  Without this the last
            # block of a section swallowed the next section's title lines.
            if line.startswith("## ") and not line.startswith("### "):
                section = line[3:].strip()
            key = None
        if line.startswith("## ") and not line.startswith("### "):
            continue
        heading = re.match(r"^#### `(.*)`", line)
        if heading:
            repeats = sum(
                1 for existing in blocks
                if existing[0] == section and existing[1] == heading.group(1)
            )
            key = (section, heading.group(1), repeats)
            blocks[key] = []
            continue
        if key is not None and line.strip():
            blocks[key].append(line.strip())
    return blocks


def to_element_basis(expression):
    """Rewrite both coordinate spellings of a poloidal bracket into one basis.

    The element routine writes the same antisymmetric product either as
    ``a_s*b_t - a_t*b_s`` or as ``xjac*(a_x*b_y - a_y*b_x)``, and it is not
    always consistent between a residual and its own tangent.  Expanding every
    ``_s``/``_t`` derivative through the chain rule,

        f_s = x_s*f_x + y_s*f_y,   f_t = x_t*f_x + y_t*f_y,
        xjac = x_s*y_t - x_t*y_s,

    maps both spellings onto the same polynomial, so a residual that vanishes
    in this basis means the two sides carry the same term.
    """

    x_s, x_t, y_s, y_t = sp.symbols("x_s x_t y_s y_t")
    replacements = {sp.Symbol("xjac"): x_s * y_t - x_t * y_s}
    for symbol in expression.free_symbols:
        name = str(symbol)
        if name in ("x_s", "x_t", "y_s", "y_t", "xjac") or name.startswith("var_"):
            continue
        if name.endswith("_s"):
            replacements[symbol] = (
                x_s * sp.Symbol(name[:-2] + "_x") + y_s * sp.Symbol(name[:-2] + "_y")
            )
        elif name.endswith("_t"):
            replacements[symbol] = (
                x_t * sp.Symbol(name[:-2] + "_x") + y_t * sp.Symbol(name[:-2] + "_y")
            )
    return sp.expand(expression.xreplace(replacements))


def to_physical(expression):
    """Render a residual with a single spelling for poloidal derivatives.

    Specialising the element map to ``x_s = y_t = 1``, ``x_t = y_s = 0`` makes
    ``f_s`` and ``f_t`` coincide with ``f_x`` and ``f_y`` and ``xjac`` with 1,
    so a residual printed this way is readable and free of the (s,t)/(R,Z)
    bookkeeping.  The verdict is always taken in the full element basis; this
    is only for display.
    """

    replacements = {sp.Symbol("xjac"): sp.Integer(1)}
    for symbol in expression.free_symbols:
        name = str(symbol)
        if name in ("xjac",) or name.startswith("var_"):
            continue
        if name.endswith("_s"):
            replacements[symbol] = sp.Symbol(name[:-2] + "_x")
        elif name.endswith("_t"):
            replacements[symbol] = sp.Symbol(name[:-2] + "_y")
    return sp.expand(expression.xreplace(replacements))


def check_alignment(source_path, generated_path, expected):
    """Verify that blank report cells account for the multiset difference.

    The reports are line aligned, so a blank generated cell marks a term the
    element routine has and the generator does not.  There must be exactly as
    many of them as the block has source-only monomials, and likewise on the
    other side.  Any surplus is a pair the matcher failed to put on one line
    although the report renders both the same way — a defect in the exporter
    that the block verdicts cannot see, because the multiset difference is
    computed from the file contents and cancels such a pair out.
    """

    source = Path(source_path).read_text(encoding="utf-8").splitlines()
    generated = Path(generated_path).read_text(encoding="utf-8").splitlines()
    section = None
    key = None
    seen = collections.Counter()
    unpaired = collections.Counter()
    for left, right in zip(source, generated):
        if left.startswith("#") or left.startswith(">"):
            if left.startswith("## ") and not left.startswith("### "):
                section = left[3:].strip()
            heading = re.match(r"^#### `(.*)`", left)
            key = None
            if heading:
                name = heading.group(1)
                key = (section, name, seen[(section, name)])
                seen[(section, name)] += 1
            continue
        if key is None or not (left.strip() or right.strip()):
            continue
        if left.strip() and not right.strip():
            unpaired[(key, "source")] += 1
        elif right.strip() and not left.strip():
            unpaired[(key, "generated")] += 1
    surplus = collections.Counter()
    for (key, side), count in unpaired.items():
        # A negative difference means two terms that do differ were still
        # placed on one line, by the coefficient or radius fallback.  That is
        # the intent; only a surplus of blanks is a failure to pair.
        extra = count - expected.get((key, side), 0)
        if extra > 0:
            surplus[(key, side)] = extra
    return surplus


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    reports = PROJECT_ROOT / "reports"
    parser.add_argument("source", nargs="?", default=reports / "model600_fortran_terms.md")
    parser.add_argument("generated", nargs="?", default=reports / "model600_generated_terms.md")
    parser.add_argument(
        "--quiet", action="store_true", help="print only the per-block verdicts",
    )
    args = parser.parse_args()

    source = read_blocks(args.source)
    generated = read_blocks(args.generated)

    differing = 0
    expected_blanks = {}
    keys = list(OrderedDict.fromkeys(list(source) + list(generated)))
    for key in keys:
        source_only = Counter(source.get(key, ())) - Counter(generated.get(key, ()))
        generated_only = Counter(generated.get(key, ())) - Counter(source.get(key, ()))
        label = "{} {} #{}".format(*key)
        expected_blanks[(key, "source")] = sum(source_only.values())
        expected_blanks[(key, "generated")] = sum(generated_only.values())
        if not source_only and not generated_only:
            print("OK   {}  ({} terms)".format(label, len(source.get(key, ()))))
            continue
        def _sum(terms):
            """Sum the parsable terms, listing the ones that are not.

            A source line the exporter could not expand is copied into the
            report verbatim; it is not a SymPy expression.
            """

            total = sp.S.Zero
            unparsed = []
            for term in terms:
                # Parse every identifier as a plain symbol: several JOREK work
                # variables (``zeta``, ``gamma``, ``beta``, ...) collide with
                # SymPy's own function names.
                locals_ = {
                    name: sp.Symbol(name)
                    for name in re.findall(r"[A-Za-z_][A-Za-z_0-9]*", term)
                    if name != "sqrt"
                }
                try:
                    total += sp.sympify(term, locals=locals_, evaluate=True)
                except Exception:
                    unparsed.append(term)
            return total, unparsed

        source_total, source_unparsed = _sum(source_only.elements())
        generated_total, generated_unparsed = _sum(generated_only.elements())
        residual = sp.expand(source_total - generated_total)
        # The verdict is taken in the element basis, where the two spellings
        # of a poloidal bracket coincide; anything left there is a genuine
        # difference.
        residual_basis = to_element_basis(residual)
        if residual != 0 and residual_basis == 0:
            # The two sides carry the same physics; the element routine simply
            # wrote a poloidal bracket in (s,t) where the generator wrote it
            # in (R,Z), or the other way round.
            print("SAME {}: {} line(s), identical up to the (s,t)/(R,Z) "
                  "spelling of a poloidal bracket".format(
                      label, sum(source_only.values()),
                  ))
            continue
        differing += 1
        print("DIFF {}: source-only={} generated-only={} residual terms (element basis)={}{}".format(
            label, sum(source_only.values()), sum(generated_only.values()),
            0 if residual_basis == 0 else len(sp.Add.make_args(residual_basis)),
            "" if not (source_unparsed or generated_unparsed)
            else " (unparsed: {})".format(
                len(source_unparsed) + len(generated_unparsed)
            ),
        ))
        if args.quiet:
            continue
        for term in sorted(source_only.elements()):
            print("   source    {}".format(term))
        for term in sorted(generated_only.elements()):
            print("   generated {}".format(term))
        # Print the residual with one spelling for the poloidal derivatives,
        # so that lines differing only in (s,t) versus (R,Z) cancel here as
        # they do in the verdict above.
        shown = to_physical(residual)
        if shown == 0:
            print("   residual: 0 in the displayed basis; the difference is "
                  "visible only in the element metric")
        else:
            print("   residual (source-generated, with f_s->f_x, f_t->f_y, "
                  "xjac->1):")
            for term in sp.Add.make_args(shown):
                print("     {}".format(term))
    print("\nBlocks with differences: {} / {}".format(differing, len(keys)))
    misaligned = check_alignment(args.source, args.generated, expected_blanks)
    if misaligned:
        print("MISALIGNED: {} report line(s) left unpaired although both "
              "reports contain the term:".format(sum(misaligned.values())))
        for (block, side), count in misaligned.items():
            print("   {} {} #{} ({} side): {}".format(*block, side, count))
    else:
        print("Alignment check: blank cells match the multiset difference in "
              "every block.")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())
