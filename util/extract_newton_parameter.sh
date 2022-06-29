#!/bin/bash
PROGRESS_BAR_WIDTH=30  # progress bar length in characters

draw_progress_bar() {
  # Arguments: current value, max value, unit of measurement (optional)
  local __value=$1
  local __max=$2
  local __unit=${3:-""}  # if unit is not supplied, do not display it

  # Calculate percentage
  if (( $__max < 1 )); then __max=1; fi  # anti zero division protection
  local __percentage=$(( 100 - ($__max*100 - $__value*100) / $__max ))

  # Rescale the bar according to the progress bar width
  local __num_bar=$(( $__percentage * $PROGRESS_BAR_WIDTH / 100 ))

  # Draw progress bar
  printf "["
  for b in $(seq 1 $__num_bar); do printf "#"; done
  for s in $(seq 1 $(( $PROGRESS_BAR_WIDTH - $__num_bar ))); do printf " "; done
  printf "] $__percentage%% ($__value / $__max $__unit)\r"
}


# loop over run folders
logfile="logfile.out"
directory=scan_folders_tstep1.0d-1
mkdir -p $directory/compare_newton
cd $directory
declare -a run_folders=(./*_tmp_newton)
length=${#run_folders[*]}
counter=0
echo "Extracting from "$length "run folders..."
shopt -s extglob
for folder in ${run_folders[@]}
do
 #echo $folder
 if [[ -e $folder/$logfile ]]; then
  dummy="$(grep 'ABORTING' $folder/$logfile)"
  if [[ $dummy == "" ]]; then 
   grep newton_i "$folder"/$logfile | tr -d " " | awk -F "t_step=|newton_i=|r_k|iter_gmres" 'BEGIN {print "# t", "  n", "  g"}{printf "%3s %3s %3s\n",$2, $3, $5}' > compare_newton/"$folder"_niterates.dat 
  cd "$folder"; rm -r !(logfile|logfile.out|intear|input); cd ..;
  fi
# else
#  echo $folder
 fi
 counter=$((counter+1))
 draw_progress_bar $counter $length
done
sleep .5
printf "\n"
printf "done."
printf "\n" 
