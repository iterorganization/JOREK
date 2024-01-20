#!/bin/bash

# The script needs to be called in a special way: '. ./util/detect_rst_type.sh'

function print_help() {
  echo "Usage: detect_rst_type.sh [options]"
  echo ""
  echo " -h        Print this help text"
  echo " -n NNNNNN  Detect restart type based on a certain time step"
  echo " -v        Print additional information"
}

if [ "$1" == "-h" ]; then
  print_help
  stop
fi

function set_rst_type() {
  export RST_TYPE="$1"
}

timestep="[0-9][0-9][0-9][0-9][0-9][0-9]"
verbose="0"
while [ $# -gt 0 ]; do
  if [ "$1" == "-n" ]; then
    timestep="$2"
    shift 2
  elif [ "$1" == "v" ]; then
    verbose="1"
    shift
  else
    echo "ERROR: Unknown option $1."
    print_help
    stop
  fi
done

num_rst=`ls -1d jorek${timestep}.rst 2>/dev/null | wc -l`
num_h5=`ls -1d jorek${timestep}.h5 2>/dev/null | wc -l`

set_rst_type "h5" # default
if [ $num_rst -gt 0 ] && [ $num_h5 -gt 0 ]; then
  echo "WARNING in detect_rst_type: Both '.rst' and '.h5' restart files exist!"
  if [ $num_rst -gt $num_h5 ]; then
    set_rst_type "rst"
  fi
elif [ $num_rst -gt 0 ]; then
  set_rst_type "rst"
elif [ $num_h5 -gt 0 ]; then
  set_rst_type "h5"
else
  echo "WARNING in detect_rst_type: No '.rst' and '.h5' restart files exist!"
fi

echo "detect_rst_type: $RST_TYPE"
