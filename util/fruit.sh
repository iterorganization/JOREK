#!/usr/bin/env bash

# 
# Compile and run the FRUIT unit tests in file(s) passed in argument
#
# Date: 2018-01-15
# Author: Daan van Vugt, TU Eindhoven
#
set -u

function usage() {
  echo ""
  echo "Usage: `basename $0` [-hk] [-t <type> --type=<type>] <file/dir>..."
  echo ""
  echo "  <file/dir> A file or directory containing FRUIT test (files)."
  echo "  -h         Print this usage information and exit"
  echo "  -k         Keep the generated test executable and source file (for running in GDB)"
  echo "  -t <type>  Type of output. Either none, junit or xml (default none)"
  echo "  -s <name>  Test only a single subroutine if specified"
  echo ""
  echo "The executable created will have a temporary file name."
}

has_setup=0
has_teardown=0
keep_executable=0
xml=""
junit=0
outfile="test"
test_only=""

while getopts ":hkt:s:" opt; do
  case $opt in
    h)
      usage
      exit 0
      ;;
    k)
      keep_executable=1
      ;;
    t)
      case $OPTARG in
        none)
          ;;
        junit)
          junit=1
          ;;
        xml)
          xml='_xml'
          ;;
        *)
          echo "Unknown output type $OPTARG" >&2
          usage
          exit 1
          ;;
      esac
      ;;
    s)
      test_only="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))


function set_outfile() {
  outfile=`mktemp test_XXX`
}

function cleanup() {
  rm -f "$outfile.mods" "$outfile.tests" "$outfile.setup" "$outfile.teardown"
  if [ "$keep_executable" -eq 0 ]; then
    rm -f "$outfile" "$outfile.f90"
  else
    echo "Test executable saved in $outfile"
  fi
}

function scanfile() {
  file=${1%/}
  if [ ! -f "$file" ]; then
    echo "Cannot open '$file'." >&2
    usage
    exit 1
  fi
  echo "Scanning $file"
  # Look for any lines named module *
  grep -Eo '^ *module [[:alnum:]_]*' "$file" | sed 's/module/use/g' >> $outfile.mods
  # Look for setup subroutine (global, can be only one)
  if grep -q 'subroutine \bsetup\b' "$file"; then
    if [ "$has_setup" -eq 1 ]; then
      echo "ERROR: 2 setup routines found, expect linking errors"
    else
      has_setup=1
    fi
  fi
  # Look for any subroutines called setup_* or *_setup
  grep -o 'subroutine setup_[^ ]*\|subroutine [^ ]*_setup' "$file" >> $outfile.setup
  # Look for any subroutines called test_* or *_test
  # might contain duplicates due to end subroutine. Remove those later
  grep -o 'subroutine test_[^ ]*' "$file" >> $outfile.tests
  # Look for any subroutines called teardown_* or *_teardo
  grep -o 'subroutine teardown_[^ ]*\|subroutine [^ ]*_teardown' "$file" >> $outfile.teardown

  # Look for teardown subroutine (global, can be only one)
  if grep -q 'subroutine \bteardown\b' "$file"; then
    has_teardown=1
  fi
}

function writetest() {
  echo "program $outfile" > $outfile.f90
  echo "use fruit" >> $outfile.f90 # use fruit_mpi if you want mpi support
  sort $outfile.mods | uniq >> $outfile.f90
  echo "implicit none" >> $outfile.f90
  echo "call init_fruit$xml" >> $outfile.f90
  if [ $has_setup -eq 1 ]; then
    echo "call setup" >> $outfile.f90
  fi
  sed -e 's/subroutine/call/' < $outfile.setup >> $outfile.f90
  if [ ! -z "$test_only" ]; then
    # Filter out the subroutine if we wish to test only one
    uniq $outfile.tests | grep -F "subroutine $test_only" | sed -e 's/subroutine \([^ ]*\)/call run_test_case(\1,"\1")/g' >> $outfile.f90
  else
    uniq $outfile.tests | sed -e 's/subroutine \([^ ]*\)/call run_test_case(\1,"\1")/g' >> $outfile.f90
  fi
  echo "call fruit_summary$xml" >> $outfile.f90
  echo "call fruit_finalize" >> $outfile.f90
  sed -e 's/subroutine/call/' < $outfile.teardown >> $outfile.f90
  if [ $has_teardown -eq 1 ]; then
    echo "call teardown" >> $outfile.f90
  fi
  echo "end program $outfile" >> $outfile.f90

  rm "$outfile.tests"
  rm "$outfile.mods"
}

function runtest() {
  make $outfile
  if [ $? -eq 0 ]; then
    if [ "$junit" -eq 1 ]; then
      ./$outfile > $outfile.log
      util/fruit2junit.sh $outfile.log
    else
      ./$outfile
    fi
  fi
}



if [ $# -lt 1 ]; then
  echo "Missing file/directory name." >&2
  usage
  exit 1
fi

trap cleanup EXIT
set_outfile

for file in `find $@ -maxdepth 1 -type f -name '*.f90' -not -name 'setup_*.f90'`; do
  scanfile $file
done
# automatically add setup_*.f90 in any of the folders mentioned
dirs=
for file in $@; do
  if [ -d "$file" ]; then
    dirs="$dirs $file"
  else
    dirs="$dirs $(dirname $file)"
  fi
done
for file in `find $dirs -maxdepth 1 -type f -name 'setup_*.f90'`; do
  scanfile $file
done

writetest
runtest
cleanup
