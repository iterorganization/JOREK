# --- General settings
jorekmodel="600"
modelsettings="with_vpar=.true. with_impurities=.true. with_TiTe=.true. with_neutrals=.true."
description="Shattered pellet injection case, JET-like plasma, model$jorekmodel ($modelsettings), n_tor=3, Neon."
mpitasks=2
binaries="jorek_model${jorekmodel}_3"
binaries_initial="jorek_model${jorekmodel}_1"
requiredfiles="input acd96_ne.dat plt96_ne.dat prc96_ne.dat ycd96_ne.dat ecd96_ne.dat prb96_ne.dat scd96_ne.dat zcd96_ne.dat ion_ne.dat"
extra_remote_files="acd96_ne.dat plt96_ne.dat prc96_ne.dat ycd96_ne.dat ecd96_ne.dat prb96_ne.dat scd96_ne.dat zcd96_ne.dat ion_ne.dat"


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel $modelsettings n_tor=1 n_plane=1 n_period=1                  || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                           || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                           || exit 1
    make cleanall                                                                    || exit 1
  fi
  ./util/config.sh model=$jorekmodel $modelsettings n_tor=3 n_plane=4 n_period=1                    || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel}                             || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3                             || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input nstep_n=50,50,100 tstep_n=0.01,0.2,5. restart=.f.|| exit 1
  $MPIRUN 1 ./jorek_model${jorekmodel}_1 < input | tee logfile_initial               || exit 1
  ${codedir}/util/setinput.sh input nstep_n=50,80,100 tstep_n=1.,5.,1. restart=.t.         || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee logfile_initial2      || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=1. nout=1          || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee logfile               || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-6                                                      || exit 1
}
