# --- General settings
jorekmodel="600"
options="with_vpar=.true. with_TiTe=.true. with_neutrals=.false. with_impurities=.false. with_refluid=.false."
description="Ballooning mode, simple X-point plasma, model$jorekmodel, n_tor=7 + FFT."
mpitasks=7
binaries="jorek_model${jorekmodel}_5"
binaries_initial="jorek_model${jorekmodel}_1"
requiredfiles="input0 input"
extra_remote_files=""


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1 $options         || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                           || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                           || exit 1
    make cleanall                                                                    || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=5 n_plane=8 n_period=3 $options           || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_5                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input0 nstep_n=10,10,10,5,5,5 tstep_n=1.d-3,1.d-2,1.d-1,1.d0,1.d1,2.d1 || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input0 | tee logfile_initial      || exit 1
  ${codedir}/util/setinput.sh input nstep_n=175 tstep_n=2.d1 restart=.t.             || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_5 < input | tee logfile_initial2      || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=5 nout=1           || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_5 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 3.e-8                                                      || exit 1
}
