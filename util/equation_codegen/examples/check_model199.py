#!/usr/bin/env python3
"""Compare all implemented model-199 equations with the Fortran source."""

from pathlib import Path
import argparse
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.fortran_source import (  # noqa: E402
    EQUATION_1_ASSIGNMENTS,
    EQUATION_2_ASSIGNMENTS,
    EQUATION_3_ASSIGNMENTS,
    EQUATION_4_ASSIGNMENTS,
    EQUATION_5_ASSIGNMENTS,
    EQUATION_6_ASSIGNMENTS,
    FortranComparisonError,
    compare_assignment_maps,
    extract_fortran_assignments,
    generated_assignment_text,
    generated_model199_transport_assignments,
    generated_partial_model199_assignments,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keep-BigR-x", dest="set_BigR_x_one", action="store_false",
        help="retain BigR_x geometry factors during comparison",
    )
    arguments = parser.parse_args()
    source_file = REPOSITORY_ROOT / "models/model199/mod_elt_matrix_fft.f90"
    assignment_names = (
        EQUATION_1_ASSIGNMENTS + EQUATION_2_ASSIGNMENTS + EQUATION_3_ASSIGNMENTS
        + EQUATION_4_ASSIGNMENTS + EQUATION_5_ASSIGNMENTS + EQUATION_6_ASSIGNMENTS
    )
    source = extract_fortran_assignments(source_file, assignment_names)
    generated = generated_assignment_text()
    generated.update(generated_partial_model199_assignments())
    generated.update(generated_model199_transport_assignments())
    try:
        compare_assignment_maps(
            source, generated, set_BigR_x_one=arguments.set_BigR_x_one
        )
    except FortranComparisonError as error:
        print(error)
        raise SystemExit(1)
    print("No differences found for model-199 equations 1–6")


if __name__ == "__main__":
    main()
