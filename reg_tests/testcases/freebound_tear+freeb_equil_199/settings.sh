# --- General settings
jorekmodel="199"
description="2/1 mode with freebnd_equil (n=0,1 free). Same growth rate 2.29d-4 if freebound_equil=.f. Slightly modified version of benchmark with CASTOR."
mpitasks=2
binaries="jorek_model${jorekmodel}_1"
binaries_initial=""
requiredfiles="input starwall-response.dat"
extra_remote_files="starwall-response.dat"


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=3 n_plane=4 n_period=1                    || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=40 tstep_n=1000.d0           || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile_initial       || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=1000.d0            || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 3.e-5                                                      || exit 1
}
