#!/bin/bash

#
# Purpose: Produces .vtk files that can be plotted, e.g., by Paraview from JOREK restart files.
#
# Date: 2011-08-22
# Author: Matthias Hoelzl, IPP Garching
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
    rm -rf $local_tmp_dir
  fi
  echo "...done."
  exit $1
}

function usage () {
  echo ""
  echo "Convert JOREK restart files into 2D/3D VTK files for visualization."
  echo ""
  echo "Usage: `basename $0` [options] binary infile [extra-files]"
  echo ""
  echo "Options:"
  echo "  -dir <dir>                  Specify a custom target directory (see remarks below)"
  echo "  -j <nthreads>               Convert using <nthreads> threads [default: 1]"
  echo "  -only <step>,<step>          Convert only listed time steps"
  echo "  -only <step>-<step>          Convert only time steps in the given range"
  echo "  -only <step>-<dstep>-<step>  Convert only time steps in the given range with given interval"
  echo "  -donly <dstep>               Equivalent to -only 0-<dstep>-99999"
  echo "  -time <time>,<time>          Selects time step roughly at <time> (JOREK-units)"
  echo "  -time <time>-<dtime>-<time>  Selects time step within given time range with given interval"
  echo "  -dtime <dtime>               Equivalent to -time 0-<dtime>-inf"
  echo "  -zip                        Compress the .vtk files using gzip"
  echo ""
  echo "Options passed to the binary via namelist input (source code for details):"
  echo "  -i_tor <i_tor>              Select one toroidal mode [default: -1] (2D VTK ONLY)"
  echo "  -i_plane <i_plane>          Select the toroidal plane [default: 1] (2D VTK ONLY)"
  echo "  -no0                        Don't include the n=0 mode for i_tor=-1 [default: off]"
  echo "  -nsub <nsub>                Number of finite element subdivisions [default: 5]"
  echo "  -[no]periodic               (Un)set periodicity for 3D plot (3D VTK ONLY)"
  echo "  -si                         Convert to VTK file with SI units (2D VTK ONLY)"
  echo "  -fluxes                     Include energy and density fluxes [default: off] (2D VTK ONLY)"
  echo "  -neo                        Include neoclassical and more terms [default: off] (2D VTK ONLY)"
  echo "  -Bfield                     Include vector of magnetic field [default: off] (2D VTK ONLY)"
  echo "  -Vfield                     Include vector of velocity field [default: off] (2D VTK ONLY)"
  echo "  -[no]psiN                   Include normalized poloidal flux or not [default: on] (2D VTK ONLY)"
  echo "  -bootstrap                  Include bootstrap current decomposition [default: off] (2D VTK ONLY)"
  echo ""
  echo "  binary                      executable (jorek2vtk, jorek2vtk_3d, jorek2_target2vtk)"
  echo "  infile                      Input file of the corresponding JOREK run"
  echo "  extra-files                 Additional files that are required for running"
  echo ""
  echo "Remarks:"
  echo "  * <i_tor> and <i_plane> may contain lists and/or ranges, e.g., -i_tor 3,5-7"
  echo "  * The default target directory is ./vtk(3d)"
  echo "  * With option -i_tor X, '_itorX' is added to the output directory name"
  echo "  * With option -no0, '_no0' is added to the output directory name"
  echo "  * With option -i_plane X, '_iplaneX' is added to the output directory name"
  echo "  * CURRENTLY MOST OPTIONS DO NOT WORK FOR jorek2_target2vtk"
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
  targetFile="jorek.$stepnum.vtk" # Target filename with same number as source
  targetFile="$targetDir/$targetFile" # Target filename with full path

  # Convert only new, selected restart files
  if ( [ ! -e $targetFile ] || [ "$file" -nt "$targetFile" ] ) \
    &&  ( [ ! -z "$select_arguments" ] || [ `is_selected $stepnum` == "yes" ] ) ; then
    rm -f jorek_restart.${RST_TYPE}
    ln -s $file jorek_restart.${RST_TYPE}
    for copyfile in $copyfiles; do
      cp $startDir/$copyfile .
    done
    $binary < $infile > ./log 2>&1
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
    mv jorek_tmp.vtk $targetFile
    if [ "$zipfiles" == "yes" ]; then
      rm -f ${targetFile}.gz
      gzip $targetFile
    fi
  fi
  
  unmark_running $ithread
}



SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`



# --- Process command line parameters
nthreads="1"
selected_steps="0-99999"
select_arguments=""
customdir=""
nsub=""
i_tor=""
i_plane=""
no0=""
periodic=""
writenml="no"
si_units=""
include_fluxes=""         # include energy and density fluxes (or not)
include_neo=""            # include neoclassical and more terms (or not)
include_magnetic_field="" # include vector of magnetic field (or not)
include_velocity_field="" # include vector of velocity field (or not)
include_electric_field="" # include vector of electric field (or not)
include_Jpol=""           # include vector of poloidal currents (or not)
include_bootstrap=""      # include bootstrap current and averaged current
include_psi_norm=".true." # include normalized flux
while [ $# -gt 1 ]; do
  if [ "$1" == "-j" ]; then
    nthreads="$2"
    shift 2
  elif [ "$1" == "-only" ]; then
    selected_steps="$2"
    #select_arguments="$1 $2"
    shift 2
  elif [ "$1" == "-donly" ]; then
    selected_steps="0-$2-99999"
    #select_arguments="$1 $2"
    shift 2
  elif ( [ "$1" == "-time" ] || [ "$1" == "-dtime" ] ) ; then
    select_arguments="$1 $2"
    shift 2
  elif [ "$1" == "-dir" ]; then
    customdir="$2"
    shift 2
  elif [ "$1" == "-zip" ]; then
    zipfiles="yes"
    shift 1
  elif [ "$1" == "-nsub" ]; then
    nsub="$2"
    shift 2
    writenml="yes"
  elif [ "$1" == "-i_tor" ]; then
    i_tor="$2"
    shift 2
    writenml="yes"
  elif [ "$1" == "-i_plane" ]; then
    i_plane="$2"
    shift 2
    writenml="yes"
  elif [ "$1" == "-no0" ]; then
    no0="yes"
    shift 1
    writenml="yes"
  elif [ "$1" == "-noperiodic" ]; then
    periodic=".false."
    shift 1
    writenml="yes"
  elif [ "$1" == "-periodic" ]; then
    periodic=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-si" ]; then
    si_units=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-fluxes" ]; then
    include_fluxes=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-neo" ]; then
    include_neo=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-Bfield" ] || [ "$1" == "-magnetic_field" ]; then
    include_magnetic_field=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-Vfield" ] || [ "$1" == "-velocity_field" ]; then
    include_velocity_field=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-Efield" ] || [ "$1" == "-electric_field" ]; then
    include_electric_field=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-Jpol" ] || [ "$1" == "-poloidal_currents" ]; then
    include_Jpol =".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-bootstrap" ]; then
    include_bootstrap=".true."
    shift 1
    writenml="yes"
  elif [ "$1" == "-nopsiN" ] || [ "$1" == "-nopsi_norm" ]; then
    include_psi_norm=".false."
    shift 1
    writenml="yes"
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

if [ $# -lt 2 ]; then
  echo "ERROR: Not enough parameters."
  usage
  exit 1
fi

binary=`readlink -f $1`
shift
infile=`readlink -f $1`
shift
sourceDir=`readlink -f .`
copyfiles=`grep _file $infile | grep -v '^ *!' | sed -e "s/^.*= *[\"']\(.*\)[\"'].*$/\1/" | grep -v 'none'`
copyfiles="$copyfiles $@"
for copyfile in $copyfiles; do
  if [ ! -f "$copyfile" ]; then
    echo "ERROR: Extra-file '$copyfile' does not exist."
    usage
    exit 1
  fi
done



function get_recursive_params () {
  i_tor_rec=$1
  i_plane_rec=$2
  params="-j $nthreads -only $selected_steps"
  if [ ! -z "$nsub"        ]; then params="$params -nsub    $nsub";        fi
  if [ ! -z "$no0"         ]; then params="$params -no0";                  fi
  if [ ! -z "$i_tor_rec"   ]; then params="$params -i_tor   $i_tor_rec";   fi
  if [ ! -z "$i_plane_rec" ]; then params="$params -i_plane $i_plane_rec"; fi
  params="$params $binary $infile $copyfiles"
  echo $params
}

# --- Handle range specifications (e.g., 3-7) and list specifications (e.g., 3,5,7)
#     for i_tor and i_plane via recursive calls
i_tor_array=( $(echo $i_tor | tr ',' '\n') )
if [ ${#i_tor_array[@]} -gt 1 ]; then
  for i in ${i_tor_array[@]}; do
    $SCRIPTDIR/convert2vtk.sh `get_recursive_params "$i" "$i_plane"` || exit 1
  done
  exit
fi
i_tor_array=( $(echo $i_tor | tr '-' '\n') )
if [ ${#i_tor_array[@]} -eq 2 ]; then
  for i in `seq ${i_tor_array[0]} ${i_tor_array[1]}`; do
    $SCRIPTDIR/convert2vtk.sh `get_recursive_params "$i" "$i_plane"` || exit 1
  done
  exit
fi
i_plane_array=( $(echo $i_plane | tr ',' '\n') )
if [ ${#i_plane_array[@]} -gt 1 ]; then
  for i in ${i_plane_array[@]}; do
    $SCRIPTDIR/convert2vtk.sh `get_recursive_params "$i_tor" "$i"` || exit 1
  done
  exit
fi
i_plane_array=( $(echo $i_plane | tr '-' '\n') )
if [ ${#i_plane_array[@]} -eq 2 ]; then
  for i in `seq ${i_plane_array[0]} ${i_plane_array[1]}`; do
    $SCRIPTDIR/convert2vtk.sh `get_recursive_params "$i_tor" "$i"` || exit 1
  done
  exit
fi



# --- Determine output directory
if [ ! -z "$customdir" ]; then
  dir="$customdir"
else
  if [[ "$binary" == *jorek2vtk_3d* ]]; then
    dir="./vtk3d"
    threeD="yes"
    target="no"
  elif [[ "$binary" == *jorek2_target2vtk* ]]; then
    dir="./vtkTarget"
    threeD="no"
    target="yes"
  else
    dir="./vtk"
    threeD="no"
    target="no"
  fi
  if [ ! -z "$i_tor" ]; then
    dir="${dir}_itor$i_tor"
    echo ""
    echo "****** i_tor=$i_tor ******"
  elif [ ! -z "$no0" ]; then
    dir="${dir}_no0"
  fi
  if [ ! -z "$i_plane" ]; then
    dir="${dir}_iplane$i_plane"
    echo ""
    echo "****** i_plane=$i_plane ******"
  elif [ -z "$i_tor" ] && [ "$threeD" == "no" ] && [ "$target" == "no" ]; then
    dir="${dir}_iplane1"
  fi
  if [ ! -z "$si_units" ]; then
    dir="${dir}_si"
  fi
fi
echo "Writing files to dir='$dir'."



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



# --- Create directory for vtk files
startDir=`pwd`
mkdir -p $dir || exit 1
targetDir=`readlink -f $dir`



# --- Create local temporary directory
local_tmp_dir0="tmp_vtk_$$"
mkdir -p $local_tmp_dir0
local_tmp_dir=`readlink -f $local_tmp_dir0` # absolute path
ERROR_STOP_FILE="$local_tmp_dir/ERROR_STOP"



# --- Create the namelist file with jorek2vtk parameters
if [ "$writenml" == "yes" ]; then
  vtk_nml="$local_tmp_dir0/vtk.nml"
  echo "&vtk_params"               > $vtk_nml
  if [ ! -z "$nsub" ]; then
    echo "  nsub    = $nsub"      >> $vtk_nml
  fi
  if [ ! -z "$i_tor" ]; then
    echo "  i_tor   = $i_tor"     >> $vtk_nml
  fi
  if [ ! -z "$i_plane" ]; then
    echo "  i_plane = $i_plane"   >> $vtk_nml
  fi
  if [ ! -z "$no0" ]; then
    echo "  without_n0_mode = .true."   >> $vtk_nml
  fi
  if [ ! -z "$periodic" ]; then
    echo "  periodic = $periodic" >> $vtk_nml
  fi
  if [ ! -z "$si_units" ]; then
    echo "  si_units = $si_units" >> $vtk_nml
  fi
  if [ ! -z "$include_fluxes" ]; then
    echo "  include_fluxes = $include_fluxes" >> $vtk_nml
  fi
  if [ ! -z "$include_neo" ]; then
    echo "  include_neo = $include_neo" >> $vtk_nml
  fi
  if [ ! -z "$include_magnetic_field" ]; then
    echo "  include_magnetic_field = $include_magnetic_field" >> $vtk_nml
  fi
  if [ ! -z "$include_velocity_field" ]; then
    echo "  include_velocity_field = $include_velocity_field" >> $vtk_nml
  fi
  if [ ! -z "$include_electric_field" ]; then
    echo "  include_electric_field = $include_electric_field" >> $vtk_nml
  fi
   if [ ! -z "$include_Jpol" ]; then
    echo "  include_Jpol = $include_Jpol" >> $vtk_nml
  fi
  if [ ! -z "$include_bootstrap" ]; then
    echo "  include_bootstrap = $include_bootstrap" >> $vtk_nml
  fi
  if [ ! -z "$include_psi_norm" ]; then
    echo "  include_psi_norm = $include_psi_norm" >> $vtk_nml
  fi
  echo "/"                        >> $vtk_nml
  copyfiles="$copyfiles $vtk_nml"
elif [ -f "vtk.nml" ]; then
  # If parameters -nsub, -i_tor, -i_plane were not provided, but
  # a vtk.nml file exists, include it automatically
  copyfiles="$copyfiles vtk.nml"
fi



# --- Prepare thread temporary folders
for i in `seq $nthreads`; do
  tmpdir[$i]="$local_tmp_dir/thread_$i"
  mkdir -p ${tmpdir[$i]}
done



. ${SCRIPTDIR}/detect_rst_type.sh
if [ "$RST_TYPE" != "h5" ] && [ "$RST_TYPE" != "rst" ]; then
  echo "ERROR: RST_TYPE not detected properly: $RST_TYPE"
  exit 0
fi



# --- Parallel file conversion
echo ""
if [ -z "$select_arguments" ]; then
  files=`ls $sourceDir/jorek?????.${RST_TYPE} 2> /dev/null`
else
  files=`${SCRIPTDIR}/select_restart_files.sh $select_arguments`
  if [ "$files" == "Files are missing" ] ; then
    echo "Files to select restart files are missing"
    exit 0
  fi
fi
for file in $files; do
  if [ -f "$ERROR_STOP_FILE" ]; then cleanup; fi
  ithread=`get_available_thread`
  if [ ! -f "$ERROR_STOP_FILE" ]; then
    mark_running $ithread
    do_convert $file $ithread &
  fi
done



sleep 1
cleanup 0
