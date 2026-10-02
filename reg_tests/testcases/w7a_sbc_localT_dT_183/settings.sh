# --- General settings
jorekmodel="183"
jorek_equilibrium_model="180"
description="Stellarator SBC (sheath BC), sbc_use_local_T=.false./.true., model$jorekmodel, n_tor=3."
mpitasks=2
binaries="jorek_model${jorekmodel}_1"
binaries_initial="jorek_model${jorek_equilibrium_model}_1"
requiredfiles="input input_init dc_W7A_unst_10kPa gvec2jorek.dat set_T_at_edge.py"
extra_remote_files="dc_W7A_unst_10kPa gvec2jorek.dat set_T_at_edge.py"

function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorek_equilibrium_model n_tor=3 n_coord_tor=21 l_pol_domm=9 \
      n_plane=4 n_period=5 n_coord_period=5 with_TiTe=.false. with_vpar=.true. || exit 1
    make $compilopt $debugoptions jorek_model${jorek_equilibrium_model} || exit 1
    mv jorek_model${jorek_equilibrium_model} jorek_model${jorek_equilibrium_model}_1 || exit 1
    make cleanall || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=3 n_coord_tor=21 l_pol_domm=9 \
    n_plane=4 n_period=5 n_coord_period=5 with_TiTe=.false. with_vpar=.true. || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel} || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1 || exit 1
}

function initial_run () {
  # 1) fresh GVEC equilibrium import, no restart
  ${codedir}/util/setinput.sh input_init restart=.f. nstep=1 tstep=1 || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorek_equilibrium_model}_1 < ./input_init | tee logfile_initial || exit 1

  # 2) inject controlled, known edge T (value-DOF, DC harmonic only) -- dev-time step,
  #    NOT required for routine CI restart_run(); only needed when regenerating/updating this test
  python3 set_T_at_edge.py jorek000001.h5 jorek_restart.h5 1.0e-3 || exit 1
}

function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep=1 tstep=1.0 nout=1 \
    time_evol_scheme='"implicit Euler"' || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile || exit 1
}

function compare_results () {
  compare_results_generic 1.e-8 || exit 1
}