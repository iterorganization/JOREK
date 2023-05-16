#!/bin/bash

#
# Purpose: Produces poincare files that can be plotted, e.g., by gnuplot from JOREK restart files.
#          Adapted from convert2vtk.sh.
#
# Date: 2019-10-25
# Author: Fabian Wieschollek, IPP Garching
#

# --- Cleanup things when the user presses Ctrl-C or the script finishes.
trap cleanup 1 2 3 6
function cleanup () {
  if [ "$1" != "0" ]; then
    echo ""
    echo "ABORTING."
    echo ""
  fi
  echo ""
  echo "Waiting for threads to finish..."
  wait
  if [ ! -z "$local_tmp_dir" ]; then
    rm -rf "$local_tmp_dir"
  fi
  echo "...done."
  exit $1
}

function usage () {
  echo ""
  echo "Convert JOREK restart files to Poincare for visualization."
  echo ""
  echo "Usage: `basename $0` [options] binary infile [extra-files]"
  echo ""
  echo "Options:"
  echo "  -dir <dir>                  Specify a custom target directory (see remarks below)"
  echo "  -j <nthreads>               Convert using <nthreads> threads [default: 1]"
  echo "  -only <step>,<step>         Convert only listed time steps"
  echo "  -only <step>-<step>         Convert only time steps in the given range"
  echo "  -only <step>-<dstep>-<step> Convert only time steps in the given range with given interval"
  echo "  -donly <dstep>              Equivalent to -only 0-<dstep>-99999"
  echo "  -time <time>,<time>         Selects time step roughly at <time> (JOREK-units)"
  echo "  -time <time>-<dtime>-<time> Selects time step within given time range with given interval"
  echo "  -dtime <dtime>              Equivalent to -time 0-<dtime>-infinity"
  echo "  -ms                         -time is given in milliseconds instead of in JOREK-units"
  echo "  -l                          Creates a file containing all selected timesteps and times,"
  echo "                              if parameter -(d)time is used (default:off)"
  echo "  -zip                        Compress the .dat files using gzip"
  echo "  -tool                       Using jorek2_poincare or poincare? [default:jorek2_poincare]"
  echo "  -input_poinc                Input file for jorek2_poincare / poincare [default:stpts / pncr]"
  echo ""
  echo "  binary                      executable (jorek2poincare)"
  echo "  infile                      Input file of the corresponding JOREK run"
  echo "  extra-files                 Additional files that are required for running"
  echo ""
}

function mark_running () {
  ithread="$1"
  touch ${tmpdir[$ithread]}/ISRUNNING
}

function unmark_running () {
  ithread="$1"
  rm -f ${tmpdir[$ithread]}/ISRUNNING
}

function is_running () {
  ithread="$1"
  if [ -f ${tmpdir[$ithread]}/ISRUNNING ]; then
    echo 'yes'
  else
    echo 'no'
  fi
}

function get_available_thread () {
  while true; do
    for i in `seq $nthreads`; do
      if [ `is_running $i` == 'no' ]; then
        echo "$i"
        return
      fi
    done
    sleep 1
  done
}

# Find out if a time step was selected for conversion by the user (-only option)
function is_selected () {
  step_number=`echo $1 | sed -e 's/^[0 ]*\(.*.\)$/\1/'` # remove leading zeros
  step_ranges=`echo $selected_steps | tr ',' ' '` # split selected_steps, e.g., 1-3,5-7 -> 1-3 5-7
  for step_range in $step_ranges; do
    step_numbers=(`echo $step_range | tr '-' ' '`) # split step_range, e.g., 1-3 -> 1 3
    if [[ ${#step_numbers[*]} -eq 1 && ${step_numbers[0]} -eq $step_number ]] || \
       [[ ${#step_numbers[*]} -eq 2 && ${step_numbers[0]} -le $step_number && ${step_numbers[1]} -ge $step_number ]] || \
       [[ ${#step_numbers[*]} -eq 3 && ${step_numbers[0]} -le $step_number && $(($step_number % ${step_numbers[1]})) -eq 0 && ${step_numbers[2]} -ge $step_number ]] ; then
      echo "yes" # the step is contained in selected_steps, so answer 'yes'
      return
    fi
  done
  echo "no"
}

function do_convert () {
  file="$1"
  ithread="$2"

  cd ${tmpdir[$ithread]}

  stepnum=${file##*/} # Remove directory from filename
  stepnum=${stepnum:5:5}
  if [ "$poinc_tool" == "poincare" ]; then
    targetFile0="$targetDir/poincare.$stepnum.h5" # Target filename with same number as source
  else
    targetFile0="$targetDir/poinc_R-Z.$stepnum.dat" # Target filename with same number as source
    targetFile1="$targetDir/poinc_rho-theta.$stepnum.dat" # Target filename with same number as source
  fi

  # Convert only new, selected restart files
  #   If -only flag is used, $select_arguments is empty and selection of steps is carried out below via 'is_selected'.
  #   If -time flag is used, selection of steps has happened before by python and every incoming step is accepted here.
  if ( [ ! -e $targetFile0 ] || [ "$file" -nt "$targetFile0" ] ) \
    &&  ( [ ! -z "$select_arguments" ] || [ `is_selected $stepnum` == "yes" ] ) ; then
    rm -f jorek_restart.${RST_TYPE}
    ln -s $file jorek_restart.${RST_TYPE}
    for copyfile in $copyfiles; do
      cp $startDir/$copyfile .
    done
    if [ "$poinc_tool" == "jorek2_poincare" ]; then
      if [ "$input_poinc" != "stpts" ]; then
        mv "$input_poinc" "stpts"
      fi
    else
      if [ "$input_poinc" != "pncr" ]; then
        mv "$input_poinc" "pncr"
      fi
    fi

    if [ "$poinc_tool" == "poincare" ]; then
      mpirun -n 1 $binary < $infile > ./log 2>&1
    else 
      $binary < $infile > ./log 2>&1
    fi

    if [ $? -ne 0 ]; then
      if [ ! -f $ERROR_STOP_FILE ]; then
        touch $ERROR_STOP_FILE
        cat ./log
        echo ""
        echo "ithread = $ithread"
        echo "file    = `basename $file`"
        echo ""
        echo "ERRORS OCCURED EXECUTING THE BINARY. SEE ABOVE."
      fi
      unmark_running $ithread
      return 1
    else
      egrep -i "warning|restart time" ./log
    fi
    if [ "$poinc_tool" == "jorek2_poincare" ]; then
      mv poinc_R-Z.dat $targetFile0
      mv poinc_rho-theta.dat $targetFile1
    else
      mv poincare.h5 $targetFile0
    fi

    if [ "$zipfiles" == "yes" ]; then
      rm -f ${targetFile0}.gz
      gzip $targetFile0
      if [ "$poinc_tool" == "jorek2_poincare" ]; then 
        rm -f ${targetFile1}.gz
        gzip $targetFile1
      fi
    fi 
  fi
  
  unmark_running $ithread
}



SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`



# --- Process command line parameters
nthreads="1"
customdir=""
selected_steps="0-99999"
input_poinc="stpts"
poinc_tool="jorek2_poincare"
select_arguments=""
while [ $# -gt 1 ]; do
  if [ "$1" == "-j" ]; then
    nthreads="$2"
    shift 2
  elif [ "$1" == "-only" ]; then
    selected_steps="$2"
    shift 2
  elif [ "$1" == "-donly" ]; then
    selected_steps="0-$2-99999"
    shift 2
  elif ( [ "$1" == "-time" ] || [ "$1" == "-dtime" ] ) ; then
    select_arguments="$select_arguments $1 $2"
    shift 2
  elif [ "$1" == "-ms" ]; then
    select_arguments="$select_arguments $1"
    shift 1
  elif [ "$1" == "-l" ]; then
    select_arguments="$select_arguments $1"
    shift 1
  elif [ "$1" == "-dir" ]; then
    customdir="$2"
    shift 2
  elif [ "$1" == "-zip" ]; then
    zipfiles="yes"
    shift 1
  elif [ "$1" == "-input_poinc" ]; then
    input_poinc="$2"
    shift 2
  elif [ "$1" == "-tool" ]; then
    poinc_tool="$2"
    shift 2
  elif [ "$1" == "-h" ] || [ "$1" = "--help" ]; then
    usage
    exit 0
  elif [ "${1:0:1}" == "-" ]; then
    echo "ERROR: Unknown option "$1"."
    usage
    exit 1
  else
    break
  fi
done



# --- Some parameter checks
if [ $# -lt 2 ]; then
  echo "ERROR: Not enough parameters."
  usage
  exit 1
fi
if [ ! -z "$select_arguments" ] && [[ "$select_arguments" != *"time"* ]]; then
  echo "WARNING: -l and -ms parameters will be ignored, if -(d)time is not set."
  select_arguments=""
fi
regexp_steps="^[0-9]{1,5}(-[0-9]{1,5}){0,2}(,[0-9]{1,5}(-[0-9]{1,5}){0,2})*$"
if [[ ! "$selected_steps" =~ $regexp_steps   ]]; then
  echo "ERROR: -(d)only-parameter given in wrong format."
  usage
  exit 1
fi

if [ $poinc_tool != "jorek2_poincare" ] && [ $poinc_tool != "poincare" ]; then
  echo "ERROR: Wrong poincare tool. Use either jorek2_poincare or poincare"
  usage
  exit 1
fi

if [ "$poinc_tool" == "poincare" ]; then
  if [ "$input_poinc" == "stpts" ]; then
    input_poinc="pncr"  
  fi
fi

binary=`readlink -f $1`
shift
infile=`readlink -f $1`
shift
sourceDir=`readlink -f .`
copyfiles=`grep _file $infile | grep -v '^ *!' | sed -e "s/^.*= *[\"']\(.*\)[\"'].*$/\1/" | grep -v 'none'`
copyfiles="$copyfiles $input_poinc $@"
for copyfile in $copyfiles; do
  if [ ! -f "$copyfile" ]; then
    echo "ERROR: Extra-file '$copyfile' does not exist."
    usage
    exit 1
  fi
done



# --- Determine output directory
if [ ! -z "$customdir" ]; then
  dir="$customdir"
else
  dir="./poincare_dir"
fi



# --- Some basic checks
if [ ! -f $binary ]; then
  echo "ERROR: $binary does not exist."
  usage
  exit 1
elif [ ! -f $infile ]; then
  echo "ERROR: $infile does not exist."
  usage
  exit 1
elif [ ! -d $sourceDir ]; then
  echo "ERROR: $sourceDir does not exist."
  usage
  exit 1
elif [ ! -d $targetDir ]; then
  echo "ERROR: $targetDir does not exist."
  usage
  exit 1
fi



# ---- Detect restart file type
. ${SCRIPTDIR}/detect_rst_type.sh
if [ "$RST_TYPE" != "h5" ] && [ "$RST_TYPE" != "rst" ]; then
  echo "ERROR: RST_TYPE not detected properly: $RST_TYPE"
  usage
  exit 1
fi



# --- Select files for conversion
if [ -z "$select_arguments" ]; then
  files=`ls $sourceDir/jorek?????.${RST_TYPE} 2> /dev/null`
else
  files=`${SCRIPTDIR}/select_restart_files.sh $select_arguments`
  if [ "${files:0:5}" == "ERROR" ] ; then
    echo "$files" | head -1
    usage
    exit 1
  fi
fi



# --- Create directory for poincare files
echo "Writing files to dir='$dir'."
startDir=`pwd`
mkdir -p $dir || exit 1
targetDir=`readlink -f $dir`



# --- Create local temporary directory
local_tmp_dir0="tmp_poinc_$$"
mkdir -p $local_tmp_dir0
local_tmp_dir=`readlink -f $local_tmp_dir0` # absolute path
ERROR_STOP_FILE="$local_tmp_dir/ERROR_STOP"



# --- Prepare thread temporary folders
for i in `seq $nthreads`; do
  tmpdir[$i]="$local_tmp_dir/thread_$i"
  mkdir -p ${tmpdir[$i]}
done

# --- Create a list of available selected files ---------------------------------
if [ $selected_steps == "0-99999" ]; then
  selected_available_files=$files
else

  file_available_restarts="available_restart_files.txt"
  file_selected_restarts="selected_restart_files.txt"
  
  rm -f $file_available_restarts $file_selected_restarts
  ls -1 $sourceDir/jorek?????.${RST_TYPE} > $file_available_restarts
  
  step_ranges=`echo $selected_steps | tr ',' ' '`
  for step_range in $step_ranges; do
    step_numbers=(`echo $step_range | tr '-' ' '`) # split step_range, e.g., 1-3 -> 1 3
    if [[ ${#step_numbers[*]} -eq 1 ]]; then
      istart=${step_numbers[0]};   istep=1;                    iend=${step_numbers[0]}
    elif [[ ${#step_numbers[*]} -eq 2 ]]; then
      istart=${step_numbers[0]};   istep=1;                    iend=${step_numbers[1]}
    elif [[ ${#step_numbers[*]} -eq 3 ]]; then
      istart=${step_numbers[0]};   istep=${step_numbers[1]};   iend=${step_numbers[2]}
    fi

    for i in `seq $istart $istep $iend`; do
      padnumber=`printf "%05d" $i`
      echo $padnumber >> $file_selected_restarts
    done
  done
 
  selected_available_files=`grep -f $file_selected_restarts $file_available_restarts`
  rm -f $file_available_restarts $file_selected_restarts
fi
# ------------------------------------------------------------------------------


# --- Parallel file conversion
echo ""
for file in $selected_available_files; do
  if [ -f "$ERROR_STOP_FILE" ]; then cleanup; fi
  ithread=`get_available_thread`
  if [ ! -f "$ERROR_STOP_FILE" ]; then
    mark_running $ithread
    do_convert $file $ithread &
  fi
done



sleep 1
cleanup 0
