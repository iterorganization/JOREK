#!/bin/bash
#
# Regression test for the equation_codegen tool.
#
# Unlike the JOREK non-regression tests under reg_tests/ (which compile
# binaries, run an MPI job and diff HDF5 restart files), this test is a pure
# Python check: it recomputes every model-600 discrepancy between the
# hand-written element routine (models/model600/mod_elt_matrix_fft.f90) and
# the symbolic linearization (util/equation_codegen/src/jorek_equations),
# and compares it against the frozen reference in
# reference/model600_discrepancies.json.  A change to either side that shifts
# which terms disagree — a real regression, not just a reordering — makes
# this test fail and print exactly which monomials are no longer accounted
# for.
#
# Follows the same pass/fail convention as reg_tests/run_test.sh: exit 0 and
# a line "Test 'equation_codegen' passed." on success, non-zero and a
# detailed failure report otherwise.  Meant to be run directly, or as its own
# step in CI — it does not go through reg_tests/run_test.sh.
#
# The only runtime dependency is the ``sympy`` package; final_test.py adds
# src/ to sys.path itself, so nothing needs to be "installed" as a package.
#
# Usage:
#   util/equation_codegen/run_test.sh              # check the current tree
#   util/equation_codegen/run_test.sh --update      # rewrite the reference
#                                                    # after an intended change
#   util/equation_codegen/run_test.sh --parallel    # check the 11 rows as
#                                                    # separate background
#                                                    # processes (one per
#                                                    # equation), for a CI
#                                                    # agent with several
#                                                    # cores.  Equivalent in
#                                                    # coverage to a plain
#                                                    # run; not combinable
#                                                    # with --update or
#                                                    # --equation.

set -u

ROWS="psi u rho vpar Ti Te T rhoimp rhon zj w "

TESTNAME="equation_codegen"
STARTDIR=$(readlink -f "$(dirname "$0")")
cd "$STARTDIR" || exit 1

PYTHON=python3

# --- Prefer the system Python if it already has sympy; this is the common
#     case and avoids touching venv/pip at all.
if ! $PYTHON -c "import sympy" >/dev/null 2>&1; then

    # --- Try the environment-modules ``sympy`` package next (e.g. Lmod's
    #     "module load sympy/1.14.0-gfbf-2025b" on the ITER Bamboo agents),
    #     which needs no network access and no venv at all.  ``module`` is a
    #     shell function normally set up by an interactive login shell, so a
    #     non-interactive script (like this one, run from a Bamboo task) may
    #     not have it yet; source Lmod's init script first if not.  Lmod's own
    #     scripts are not written for ``set -u``, so relax it around this.
    set +u
    # A non-interactive, non-login shell (as Bamboo runs this) does not
    # source /etc/profile.d, so MODULEPATH itself may be unset even before
    # getting to the ``module`` function below; set it up the same way a
    # login shell would.
    if [ -z "$MODULEPATH" ]; then
        for f in /etc/profile.d/*.sh; do
            [ -r "$f" ] && . "$f" >/dev/null 2>&1
        done
    fi
    if ! type module >/dev/null 2>&1; then
        for init in \
            "${LMOD_PKG:-}/init/bash" \
            /usr/share/lmod/lmod/init/bash \
            /etc/profile.d/lmod.sh \
            /usr/share/Modules/init/bash
        do
            [ -n "$init" ] && [ -r "$init" ] && . "$init" && break
        done
    fi
    if type module >/dev/null 2>&1; then
        module load sympy/1.14.0-gfbf-2025b 2>/dev/null
        if ! $PYTHON -c "import sympy" >/dev/null 2>&1; then
            # Fall back to whatever sympy module this host happens to have,
            # in case the pinned version above is not the one available here.
            sympy_module=$(module -t avail sympy 2>&1 | grep -m1 '^sympy/')
            [ -n "$sympy_module" ] && module load "$sympy_module" 2>/dev/null
        fi
    fi
    set -u
fi

if ! $PYTHON -c "import sympy" >/dev/null 2>&1; then

    # --- Otherwise, set up a virtual environment the first time this runs.
    #     --system-site-packages lets it fall back to whatever the system
    #     Python already provides (sympy included, if installed there).
    if [ ! -x ".venv/bin/python" ]; then
        echo "Setting up .venv for ${TESTNAME} ..."
        $PYTHON -m venv --system-site-packages .venv || exit 1

        # Some Python installs (e.g. RHEL's python3 without python3-pip) create
        # a venv with no pip inside it.  Bootstrap it explicitly rather than
        # assume it is there.
        if ! .venv/bin/python -m pip --version >/dev/null 2>&1; then
            .venv/bin/python -m ensurepip --upgrade >/dev/null 2>&1
        fi

        if .venv/bin/python -m pip --version >/dev/null 2>&1; then
            # A plain package install, not an editable/local one: this only
            # needs sympy on sys.path, so there is nothing to build.
            .venv/bin/python -m pip install --quiet 'sympy>=1.12,<2' || exit 1
        elif ! .venv/bin/python -c "import sympy" >/dev/null 2>&1; then
            echo "ERROR: no working pip in .venv and no system sympy available;"
            echo "       install sympy for ${PYTHON} manually and re-run."
            exit 1
        fi
    fi

    PYTHON=.venv/bin/python
fi

# --- Pull --parallel out of the arguments; everything else passes through
#     to final_test.py as before.
parallel="no"
args=()
for arg in "$@"; do
    if [ "$arg" = "--parallel" ]; then
        parallel="yes"
    else
        args+=("$arg")
    fi
done

# --- The weak-form documentation is generated from model600.py; fail early if
#     it has not been re-rendered after a change to the equations.
if ! $PYTHON examples/model600_docs.py check; then
    echo "Test '${TESTNAME}' FAILED."
    exit 1
fi

if [ "$parallel" = "yes" ]; then
    for arg in "${args[@]:-}"; do
        case "$arg" in
            --update)
                echo "ERROR: --parallel cannot be combined with --update" \
                     "(concurrent writers would race on the reference file)."
                exit 1
                ;;
            --equation)
                echo "ERROR: --parallel already checks every row;" \
                     "drop --equation."
                exit 1
                ;;
        esac
    done
    # --- Check every row as its own background process.  Each row's blocks
    #     are independent of every other row's, so this finds exactly what a
    #     single full run would; it is just faster on a multi-core agent.
    #
    #     Each row gets one CPU: at most NPROC rows run at a time (NPROC
    #     from the CPU affinity mask this process actually has, which on
    #     Linux already reflects a cgroup/container CPU quota, not just the
    #     host's total core count), and when ``taskset`` is available each
    #     one is pinned to its own core so it cannot spill onto a neighbour's.
    #     SymPy itself is single-threaded, but a couple of underlying C
    #     libraries thread on their own unless told not to, so that is
    #     disabled too.
    NPROC=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
           NUMEXPR_NUM_THREADS=1

    logdir=$(mktemp -d)
    trap 'rm -rf "$logdir"' EXIT
    cpu=0
    running=0
    for row in $ROWS; do
        (
            if command -v taskset >/dev/null 2>&1; then
                taskset -c "$((cpu % NPROC))" \
                    $PYTHON final_test.py --equation "$row" "${args[@]}" \
                    > "$logdir/$row.log" 2>&1
            else
                $PYTHON final_test.py --equation "$row" "${args[@]}" \
                    > "$logdir/$row.log" 2>&1
            fi
            echo $? > "$logdir/$row.status"
        ) &
        cpu=$((cpu + 1))
        running=$((running + 1))
        if [ "$running" -ge "$NPROC" ]; then
            wait -n
            running=$((running - 1))
        fi
    done
    wait

    status=0
    for row in $ROWS; do
        rc=$(cat "$logdir/$row.status" 2>/dev/null || echo 1)
        if [ "$rc" != "0" ]; then
            status=1
            echo "--- ${row}: FAILED ---"
            cat "$logdir/$row.log"
        else
            tail -n 1 "$logdir/$row.log"
        fi
    done
else
    # --- Run the regression check.  final_test.py prints, on failure, exactly
    #     which monomials each disagreeing block is missing or has gained.
    $PYTHON final_test.py "${args[@]}"
    status=$?
fi

if [ $status -eq 0 ]; then
    echo "Test '${TESTNAME}' passed."
else
    echo "Test '${TESTNAME}' FAILED."
fi

exit $status
