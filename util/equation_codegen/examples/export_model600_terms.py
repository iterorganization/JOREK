#!/usr/bin/env python3
"""Export aligned model-600 source/generated terms for a visual diff."""

from pathlib import Path
import argparse
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.model600_markdown import export_model600_markdown  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--equation",
        choices=("psi", "u", "zj", "w", "rho", "vpar", "rhoimp", "rhon", "Ti", "Te", "T"),
        action="append",
        help="export only this equation; repeat the option for multiple equations",
    )
    parser.add_argument(
        "--include-neo",
        action="store_true",
        help="include NEO terms (omitted by default for the base-equation comparison)",
    )
    args = parser.parse_args()
    reports = PROJECT_ROOT / "reports"
    source, generated = export_model600_markdown(
        REPOSITORY_ROOT / "models/model600/mod_elt_matrix_fft.f90",
        reports / "model600_fortran_terms.md",
        reports / "model600_generated_terms.md",
        equations=args.equation,
        include_neo=args.include_neo,
    )
    print("Wrote {}".format(source))
    print("Wrote {}".format(generated))
    if args.equation:
        print("Selected equations: {}".format(", ".join(args.equation)))
    print("Compare with: meld {} {}".format(source, generated))


if __name__ == "__main__":
    main()
