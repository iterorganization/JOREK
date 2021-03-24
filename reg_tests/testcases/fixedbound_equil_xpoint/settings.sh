# --- General settings
jorekmodel="303"
description="Fixed boundary equilibrium in X-point plasma, model$jorekmodel, n_tor=1."
mpitasks=1
binaries="jorek_model${jorekmodel}_1"
binaries_initial=""
requiredfiles="input"
extra_remote_files=""


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1                    || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  dummy_initial_run
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=0                            || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-8                                                      || exit 1
}
