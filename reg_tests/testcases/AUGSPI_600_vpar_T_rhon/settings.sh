# --- General settings
jorekmodel="600"
options="with_vpar=.true. with_TiTe=.false. with_neutrals=.true. with_impurities=.false. with_refluid=.false."
description="Deuterium SPI driven 2/1 mode in AUG with model600"
mpitasks=2
binaries="jorek_model${jorekmodel}_3"
binaries_initial="jorek_model${jorekmodel}_1"
requiredfiles="input rho.dat t.dat ffp.dat dperp.dat zkperp.dat acd89_ar.dat ccd89_ar.dat pls89_ar.dat plt89_ar.dat prb89_ar.dat prc89_ar.dat qcd89_ar.dat scd89_ar.dat xcd89_ar.dat"
extra_remote_files="rho.dat t.dat ffp.dat dperp.dat zkperp.dat acd89_ar.dat ccd89_ar.dat pls89_ar.dat plt89_ar.dat prb89_ar.dat prc89_ar.dat qcd89_ar.dat scd89_ar.dat xcd89_ar.dat"


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1 $options         || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                           || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                           || exit 1
    make cleanall                                                                    || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=3 n_plane=4 n_period=1 $options           || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input nstep_n=20,20,20,20,20 tstep_n=1.d-3,1.d-2,1.d-1,1.d0,1.d1 || exit 1
  $MPIRUN 1 ./jorek_model${jorekmodel}_1 < input | tee logfile_initial                         || exit 1
  ${codedir}/util/setinput.sh input nstep_n=20,80,400 tstep_n=1.d1,4.d0,2.d0 restart=.t.       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee logfile_initial2                || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=2.d0 nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee logfile              || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 3.e-8                                                     || exit 1
}
