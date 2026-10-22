#!/usr/bin/env python3
"""Recompute every model-600 discrepancy and compare it with the frozen reference.

The reference records, for every ``rhs_ij``/``amat`` assignment of the element
routine, which monomials the element routine has and the generated
linearization does not, and the other way round.  Those two lists are the
complete statement of what disagrees; everything else in the reports is
identical on both sides.

Run it with no arguments to check the current tree:

    ./final_test.py

It exports the reports into a temporary directory, so it never touches
``reports/``.  It exits 0 when the discrepancies match the reference and 1
when they do not, printing what moved.  ``--update`` rewrites the reference
after an intended change; annotations of blocks that still exist are carried
over.
"""

import argparse
import hashlib
import json
import re
import sys
import tempfile
from collections import Counter, OrderedDict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "examples"))

from diff_model600_reports import read_blocks  # noqa: E402
from jorek_equations.model600_markdown import export_model600_markdown  # noqa: E402

REFERENCE = PROJECT_ROOT / "reference" / "model600_discrepancies.json"
SOURCE = REPOSITORY_ROOT / "models" / "model600" / "mod_elt_matrix_fft.f90"
FORMAT_VERSION = 1


def block_name(key):
    """Render a block key as the string used in the reference file."""

    section, lhs, repeat = key
    return "{} | {} | #{}".format(section, lhs, repeat)


def _row_of(name):
    """The equation row a block name belongs to, e.g. ``ti`` in amat(var_ti,..)."""

    match = re.search(r"\(var_([a-z]+)", name.split("|")[1])
    return match.group(1) if match else ""


def measure(source_path=SOURCE, rows=None):
    """Export the reports and return the discrepancies of every block."""

    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        source_report, generated_report = export_model600_markdown(
            source_path, folder / "source.md", folder / "generated.md",
            equations=rows,
        )
        source = read_blocks(source_report)
        generated = read_blocks(generated_report)
    blocks = OrderedDict()
    for key in list(dict.fromkeys(list(source) + list(generated))):
        source_only = Counter(source.get(key, ())) - Counter(generated.get(key, ()))
        generated_only = Counter(generated.get(key, ())) - Counter(source.get(key, ()))
        if not source_only and not generated_only:
            continue
        blocks[block_name(key)] = {
            "source_only": sorted(source_only.elements()),
            "generated_only": sorted(generated_only.elements()),
        }
    return blocks


def build(blocks, source_path=SOURCE, annotations=None):
    """Wrap the measured blocks with the metadata stored in the reference."""

    annotations = annotations or {}
    digest = hashlib.sha256(Path(source_path).read_bytes()).hexdigest()
    document = OrderedDict()
    document["format_version"] = FORMAT_VERSION
    document["source"] = str(Path(source_path).relative_to(REPOSITORY_ROOT))
    document["source_sha256"] = digest
    document["block_count"] = len(blocks)
    document["source_only_lines"] = sum(
        len(item["source_only"]) for item in blocks.values()
    )
    document["generated_only_lines"] = sum(
        len(item["generated_only"]) for item in blocks.values()
    )
    document["blocks"] = OrderedDict(
        (name, OrderedDict((
            ("findings", annotations.get(name, [])),
            ("source_only", item["source_only"]),
            ("generated_only", item["generated_only"]),
        )))
        for name, item in sorted(blocks.items())
    )
    return document


def compare(reference, blocks):
    """Return a list of human-readable differences, empty when they agree."""

    stored = reference.get("blocks", {})
    problems = []
    for name in sorted(set(stored) | set(blocks)):
        want = stored.get(name)
        have = blocks.get(name)
        if want is None:
            problems.append(
                "NEW BLOCK      {}: {} source-only, {} generated-only".format(
                    name, len(have["source_only"]), len(have["generated_only"]),
                )
            )
            problems.extend("      + " + t for t in have["source_only"])
            problems.extend("      + " + t for t in have["generated_only"])
            continue
        if have is None:
            problems.append(
                "BLOCK NOW AGREES {} (the reference expects {} source-only and "
                "{} generated-only)".format(
                    name, len(want["source_only"]), len(want["generated_only"]),
                )
            )
            continue
        for side in ("source_only", "generated_only"):
            gone = Counter(want[side]) - Counter(have[side])
            extra = Counter(have[side]) - Counter(want[side])
            if not gone and not extra:
                continue
            problems.append("{} [{}] findings {}".format(
                name, side, want.get("findings") or "unannotated",
            ))
            problems.extend("      - " + t for t in sorted(gone.elements()))
            problems.extend("      + " + t for t in sorted(extra.elements()))
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--reference", type=Path, default=REFERENCE,
        help="reference file to compare against (default: %(default)s)",
    )
    parser.add_argument(
        "--source", type=Path, default=SOURCE,
        help="element routine to read (default: the model-600 one)",
    )
    parser.add_argument(
        "--equation", action="append", dest="rows", metavar="ROW",
        help="check only this row; repeat for several.  A partial run compares "
             "only the blocks of the rows it exported, which is meant for "
             "iterating, not for signing off a change",
    )
    parser.add_argument(
        "--update", action="store_true",
        help="rewrite the reference from the current tree",
    )
    arguments = parser.parse_args()

    if arguments.rows and arguments.update:
        parser.error("--update needs a full run; drop --equation")
    print("Exporting model-600 terms and collecting discrepancies ...")
    blocks = measure(arguments.source, arguments.rows)
    source_only = sum(len(item["source_only"]) for item in blocks.values())
    generated_only = sum(len(item["generated_only"]) for item in blocks.values())
    print("  {} block(s) disagree: {} source-only and {} generated-only "
          "monomials".format(len(blocks), source_only, generated_only))

    if arguments.update:
        previous = {}
        if arguments.reference.exists():
            previous = json.loads(arguments.reference.read_text(encoding="utf-8"))
        annotations = {
            name: item.get("findings", [])
            for name, item in previous.get("blocks", {}).items()
        }
        document = build(blocks, arguments.source, annotations)
        arguments.reference.parent.mkdir(parents=True, exist_ok=True)
        arguments.reference.write_text(
            json.dumps(document, indent=1) + "\n", encoding="utf-8",
        )
        missing = [n for n in document["blocks"] if not document["blocks"][n]["findings"]]
        print("Wrote {}".format(arguments.reference))
        if missing:
            print("  {} block(s) carry no finding annotation:".format(len(missing)))
            for name in missing:
                print("    {}".format(name))
        return 0

    if not arguments.reference.exists():
        print("No reference at {}; create one with --update.".format(
            arguments.reference))
        return 1
    reference = json.loads(arguments.reference.read_text(encoding="utf-8"))

    digest = hashlib.sha256(Path(arguments.source).read_bytes()).hexdigest()
    if digest != reference.get("source_sha256"):
        print("NOTE: {} has changed since the reference was written.".format(
            reference.get("source")))
        print("      reference sha256 {}".format(reference.get("source_sha256")))
        print("      current   sha256 {}".format(digest))

    if arguments.rows:
        wanted = {row.lower() for row in arguments.rows}
        reference = dict(reference)
        reference["blocks"] = {
            name: item for name, item in reference["blocks"].items()
            if _row_of(name) in wanted
        }
        print("  comparing the {} recorded block(s) of row(s) {}".format(
            len(reference["blocks"]), ", ".join(sorted(wanted))))
    problems = compare(reference, blocks)
    if problems:
        print("\nFINAL TEST: FAILED — {} block(s) differ from the "
              "reference".format(sum(1 for p in problems if not p.startswith("   "))))
        for line in problems:
            print("  " + line)
        print("\nA '-' line is a discrepancy the reference expects but the run no "
              "longer produces; a '+' line is a new one.")
        return 1
    print("\nFINAL TEST: PASSED — every block reproduces the recorded "
          "discrepancies exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
