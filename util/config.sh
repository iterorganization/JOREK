#!/bin/bash

#
# Purpose: Set the model in the makefile and/or certain parameters in the respective
#   settings files.
#
# Date: 2011-2021
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Purpose: Modify or print physics model in Makefile.inc and further"
  echo " parameters in the corresponding mod_settings and mod_model_settings files."
  echo ""
  echo "Usage: `basename $0` [<key1>=<value1> [...]]   Modify model and/or parameters"
  echo "       `basename $0` -p <key>                  Print the value for <key> and exit"
  echo ""
  echo "Examples:"
  echo "  `basename $0` model=302 n_tor=3 n_period=8 n_plane=4"
  echo "  `basename $0` -p model"
  echo ""
}

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage && exit
elif [ -z "Makefile.inc" ]; then
  echo "Could not find a ./Makefile.inc. Are you in the JOREK trunk?" >&2
  exit 1
fi

function key() {
  echo $1 | sed -e 's/=.*$//' 
}

function val() {
  echo $1 | sed -e 's/^.*=//' 
}

function setmodel() {
  model=$1
  # --- Some checks
  if [ -e "models/model$model" ]; then
    model="model$model"
  elif [ ! ${#model} -eq 8 ] || [[ ! ${model:5:3} =~ ^[0-9]+$ ]]; then
    echo "ERROR: Illegal model specified: '$model'." >&2
    exit 1
  fi
  # --- Set model in makefile configuration files
  sed -i -e "s/\(^ *MODEL *= *\)[^ ]*\(.*$\)/\1$model\2/" Makefile.inc
  # --- Clean up .d/.o/.mod files because of the model change
  make cleanall
}

function getmodel() {
  egrep "MODEL *= *model[0-9]*" Makefile.inc | sed -e "s/^ *MODEL *= *\(model[0-9]*\).*$/\1/"
}

function setparam() {
  key=$1
  val=$2
  matches1=`grep -c ":: *$key" $paramfile1`
  matches2=`grep -c ":: *$key" $paramfile2`
  if [ "$matches1" -ne 1 ] && [ "$matches2" -ne 1 ]; then
    echo "ERROR: Could not set parameter $key." >&2
    exit 1
  else
    for paramfile in $paramfiles; do
      sed -i -e "s/\(^.*:: *$key *= *\)[^ !\t]*\(.*$\)/\1$val\2/" $paramfile
    done
  fi
}

function getparam() {
  key=$1
  for i in $paramfiles; do
    grep -i ":: *$key[ =]" $i | sed -e "s/^.*:: *$key *= *\([^ !\t]*\).*$/\1/i"
  done
}

function whichparams() {
  for i in $@; do
    grep '#SETTINGS#' $@ | sed -e 's/^.*#//'
  done
}

function print_info() {
  echo ""
  echo "=============================="
  
  params1=`whichparams $paramfile1`
  params2=`whichparams $paramfile2`
  
  echo "  `getmodel`"
  echo "------------------------------"
  for param in $params1; do
    printf "  %-16s = %s\n" "$param" "`getparam $param`"
  done
  if [ "$params2" != "" ]; then
    echo "------------------------------"
    for param in $params2; do
      printf "  %-16s = %s\n" "$param" "`getparam $param`"
    done
  fi
  echo "=============================="
  echo ""
}

function check_param_files() {
  for i in $paramfiles; do
    if [ ! -f $i ]; then
      echo "ERROR: File '$i' does not exist." >&2
      exit 1
    fi
  done
}

# --- Determine the model
model=`getmodel`
paramfile1="models/mod_settings.f90"
paramfile2="models/$model/mod_model_settings.f90"
paramfiles="$paramfile1 $paramfile2"
check_param_files

# --- If argument -p is given, just print the requested parameter value and exit
if [ "$1" == "-p" ]; then
  if [ "$2" == "model" ]; then
    echo `getmodel | sed -e 's/model//'`
  else
    value=`getparam $2`
    if [ -z "$value" ]; then
      echo "Could not find parameter '$2'."
      exit
    fi
    echo "$value"
  fi
  exit
fi

# --- First set the model (if it is specified as a command line argument)
for arg in $@; do
  if [ `key $arg` == "model" ]; then
    setmodel `val $arg`
  fi
done
model=`getmodel`
paramfile1="models/mod_settings.f90"
paramfile2="models/$model/mod_model_settings.f90"
paramfiles="$paramfile1 $paramfile2"
check_param_files

# --- Set the parameters
for arg in $@; do
  if [ `key $arg` != "model" ]; then
    setparam `key $arg` `val $arg`
  fi
done

# --- Print the configuration
model=`getmodel`
print_info $file
echo "('`basename $0` -h' for help)"
