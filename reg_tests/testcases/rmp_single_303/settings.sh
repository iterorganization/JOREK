# --- General settings
jorekmodel="303"
description="RMP with flows, n_tor=3"
mpitasks=2
binaries="jorek_model${jorekmodel}_3"
binaries_initial="jorek_model${jorekmodel}_1"
requiredfiles="input RMP_psi_cos_JET_N3.txt RMP_psi_sin_JET_N3.txt RMP_start_time.dat"
extra_remote_files="RMP_psi_cos_JET_N3.txt RMP_psi_sin_JET_N3.txt RMP_start_time.dat"


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=12                 || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                           || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                           || exit 1
    make cleanall                                                                    || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=3 n_plane=4 n_period=3                    || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ./jorek_model${jorekmodel}_1 < input | tee logfile_initial                         || exit 1
  ${codedir}/util/setinput.sh input restart=.t. tstep_n=5. nstep_n=222 iter_precon=0 gmres_4=1.d4 RMP_on=.t. RMP_psi_cos_file=\'RMP_psi_cos_JET_N3.txt\' RMP_psi_sin_file=\'RMP_psi_sin_JET_N3.txt\' nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input  | tee logfile_initial2     || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. tstep_n=5. nstep_n=1 iter_precon=0 gmres_4=1.d4 RMP_on=.t. RMP_psi_cos_file=\'RMP_psi_cos_JET_N3.txt\' RMP_psi_sin_file=\'RMP_psi_sin_JET_N3.txt\' nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-7                                                      || exit 1
}
