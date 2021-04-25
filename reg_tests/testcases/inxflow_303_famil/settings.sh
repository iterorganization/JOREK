# --- General settings
jorekmodel="303"
description="Ballooning mode, simple X-point plasma, model$jorekmodel, n_tor=7 + FFT + STRUMPACK."
mpitasks=8
binaries="jorek_model${jorekmodel}_7"
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
    ./util/config.sh model=$jorekmodel n_tor=7 n_plane=16 n_period=2                 || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                           || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_7                           || exit 1
    make cleanall                                                                    || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=7 n_plane=16 n_period=2                   || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_7                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=5,5,5,5,5 tstep_n=0.1,1,5,10,20          || exit 1
  $MPIRUN 1 ./jorek_model${jorekmodel}_1 < input | tee logfile_initial               || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=70 tstep_n=20                || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_7 < input | tee logfile_initial2      || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=20                 || exit 1
  ${codedir}/util/setinput.sh input centralize_harm_mat=.f.                          || exit 1
  ${codedir}/util/setinput.sh input autodistribute_modes=.f.                         || exit 1
  ${codedir}/util/setinput.sh input autodistribute_ranks=.f.                         || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_7 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-8                                                     || exit 1
}
