# --- General settings
jorekmodel="303"
description="Tearing mode, circular plasma, model$jorekmodel, n_tor=3 + FFT."
mpitasks=2
binaries="jorek_model${jorekmodel}_3"
binaries_initial="jorek_model${jorekmodel}_1"
requiredfiles="input"
extra_remote_files=""


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1                  || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                           || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                           || exit 1
    make cleanall                                                                    || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=3 n_plane=4 n_period=1                    || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input nstep_n=10,10,10, tstep_n=1.,100.,3000.          || exit 1
  $MPIRUN 1 ./jorek_model${jorekmodel}_1 < input | tee logfile_initial               || exit 1
  ${codedir}/util/setinput.sh input nstep_n=30 tstep_n=3000. restart=.t.             || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee logfile_initial2      || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=3000. nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-12                                                     || exit 1
}
