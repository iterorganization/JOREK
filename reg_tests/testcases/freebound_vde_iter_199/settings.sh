# --- General settings
jorekmodel="199"
description="VDE test case for an ITER plasma with simplified wall geometry (JOREK-STARWALL, model199)."
mpitasks=1
binaries="jorek_model${jorekmodel}_1"
binaries_initial=""
requiredfiles="input starwall-response.dat coil_field.dat"
extra_remote_files="starwall-response.dat coil_field.dat"


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1                    || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=5,5,10,10,90,10,170 tstep_n=1.d0,1.d1,1.d2,3.d2,1.d3,3.d2,1.d2 || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile_initial       || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=1.d2               || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-5                                                      || exit 1
}
