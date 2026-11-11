#!/usr/bin/env python3
"""Keep the model-600 weak-form documentation in step with model600.py.

    model600_docs.py render   # rewrite the generated regions of the page
    model600_docs.py check    # change nothing; exit 1 if the page is stale

The equations are written in ``src/jorek_equations/model600.py``; the page
``docs/physics/base_fluid_models/RMHD/weak_form.md`` is generated from them.
"""

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from jorek_equations.model600_docs import render_page  # noqa: E402

PAGE = REPOSITORY_ROOT / "docs" / "physics" / "base_fluid_models" / "RMHD" / "weak_form.md"
SOURCE = PROJECT_ROOT / "src" / "jorek_equations" / "model600.py"
REFERENCE = PROJECT_ROOT / "reference" / "model600_discrepancies.json"


def rendered():
    text = PAGE.read_text(encoding="utf-8")
    return text, render_page(
        text,
        SOURCE.read_text(encoding="utf-8"),
        REFERENCE.read_text(encoding="utf-8"),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("command", choices=("render", "check"))
    arguments = parser.parse_args()
    name = PAGE.relative_to(REPOSITORY_ROOT)
    text, new = rendered()
    if arguments.command == "render":
        if new != text:
            PAGE.write_text(new, encoding="utf-8")
            print("Rewrote the generated regions of {}".format(name))
        else:
            print("{} is up to date".format(name))
        return 0
    if new != text:
        print("STALE: {} does not match model600.py; run 'python3 examples/model600_docs.py "
              "render'".format(name))
        return 1
    print("{} is in step with model600.py".format(name))
    return 0


if __name__ == "__main__":
    sys.exit(main())
