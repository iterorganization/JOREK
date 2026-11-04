# --- General settings
jorekmodel="183"
jorek_equilibrium_model="180"
description="Ballooning mode, classical stellarator, model$jorekmodel, n_tor=5."
mpitasks=9
binaries="jorek_model${jorekmodel}_1"
binaries_initial="jorek_model${jorek_equilibrium_model}_1"
requiredfiles="input input_init dc_W7A_unst_10kPa gvec2jorek.dat"
extra_remote_files="dc_W7A_unst_10kPa gvec2jorek.dat"


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorek_equilibrium_model n_tor=5 n_coord_tor=21 l_pol_domm=9 n_plane=8 n_period=5 n_coord_period=5 with_TiTe=.true. with_vpar=.false. || exit 1
    make $compilopt $debugoptions jorek_model${jorek_equilibrium_model}                                                        || exit 1
    mv jorek_model${jorek_equilibrium_model} jorek_model${jorek_equilibrium_model}_1                                           || exit 1
    make cleanall                                                                                                              || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=5 n_coord_tor=21 l_pol_domm=9 n_plane=8 n_period=5 n_coord_period=5 with_TiTe=.true. with_vpar=.false.              || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                                                                       || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                                                                       || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input_init restart=.f. nstep=1 tstep=1                                          || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorek_equilibrium_model}_1 < ./input_init | tee logfile_initial              || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep=30 tstep=1.0 time_evol_scheme='"implicit Euler"'          || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < ./input | tee logfile_prerun                                 || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep=80 tstep=1.0 time_evol_scheme='"Crank-Nicholson"'         || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < ./input | tee logfile_linear                                  || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep=30 tstep=0.3 time_evol_scheme='"Crank-Nicholson"'         || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < ./input | tee logfile_nonlinear                                 || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep=1 tstep=0.3 nout=1  time_evol_scheme='"Crank-Nicholson"'  || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile                                           || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-8                                                     || exit 1
}
