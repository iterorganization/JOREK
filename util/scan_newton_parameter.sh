#!/bin/bash
start=$(date +%s.%N)
PROGRESS_BAR_WIDTH=50  # progress bar length in characters

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

current=0
logfile="logfile.out"
start_ar=(0.1 0.5 1.0 2.0) #($(seq 0.55 0.1 0.55))
alpha_ar=($(seq 0.2 0.2 2.0))
gamma_ar=($(seq 0.1 0.1 1.0))
eps_0_ar=(1.d-3) #(1.d-6 1.d-5 1.d-4 1.d-3)
length=$[${#start_ar[@]}*${#alpha_ar[@]}*${#gamma_ar[@]}*${#eps_0_ar[@]}]

for newton_start in "${start_ar[@]}"
do
 for newton_alpha in "${alpha_ar[@]}"
 do
  for newton_gamma in  "${gamma_ar[@]}"
  do
   for newton_eps_0 in "${eps_0_ar[@]}"
   do
    directory=scan_folders_tstep1.0d-1/"$newton_start"_"$newton_alpha"_"$newton_gamma"_"$newton_eps_0"_tmp_newton
    if [[ ! -e $directory/$logfile ]]; then
     mkdir -p $directory;
     cp -a -u run/. $directory;
     cd "$directory";
     ./util/setinput.sh input tstep=.1 nstep=1 newton_start=$newton_start newton_alpha=$newton_alpha newton_gamma=$newton_gamma  newton_eps_0=$newton_eps_0 newton_eps_gmres=1.d-7 restart=.t.| grep abc;
     # last line to suppress output of ./util/setinput.sh
     sbatch jobscript | grep abc;
     cd ../..;
    fi
    #else
    #echo "$directory already exists"
    current=$((current+1))
    draw_progress_bar $current $length
    # sleep 5
   done
  done
 done
done

dur=$(echo "$(date +%s.%N) - $start"|bc)
printf "\nExecution time %.6f\n" $dur
nojR=$(squeue -u jrein|grep -P -o -i " R "|wc -l)
nojPD=$(squeue -u jrein|grep -P -o -i " PD "|wc -l)
printf "$nojR jobs running, $nojPD jobs pending\n"




#squeue -u jrein
