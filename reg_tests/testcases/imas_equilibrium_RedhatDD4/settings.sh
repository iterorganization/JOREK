# --- General settings
jorekmodel="600"
options="with_vpar=.true. with_TiTe=.false. with_neutrals=.false. with_impurities=.false. with_refluid=.false. with_rho=.true."
description="Exports restart file to IMAS, then reads it from IMAS to create a JOREK input file and compares it with a reference."
mpitasks=1
binaries="jorek2_IDS" 
binaries_initial="jorek2_IDS" 
python_scripts="./communication/IMAS/imas2jorek.py" 
requiredfiles="input0 imas.nml jorek_namelist_ref"
extra_remote_files=""

# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1 $options           || exit 1
  make $compilopt $debugoptions jorek2_IDS                                           || exit 1
}

# --- Initial run only required when preparing or updating the test case
function initial_run () {
  echo "Initial run is not necessary for this case"                                  || exit 1
}

# --- Carry out the test case
function restart_run () {
   database="imas_regtest_db"
   run_number=$(shuf -i 0-100   -n 1)
   if [ -d ~/public/imasdb/${database}/4/111111/${run_number}/ ]; then
       echo "Directory exists. Delete IDS"
       rm -rf ~/public/imasdb/${database}/4/111111/${run_number}
   fi

   echo $run_number > "IMAS_RUN_OUT"
   echo $database > "IMAS_DB_OUT"
   cp jorek_restart.h5 jorek00000.h5
   sed -i "s/run_number_replace/$run_number/" imas.nml
  ./jorek2_IDS < input0                                                              || exit 1
  python imas2jorek.py -d ${database} -p 111111 -r $run_number -tk inxflow     || exit 1
  if [ -d ~/public/imasdb/${database}/4/111111/${run_number}/ ]; then
      rm -rf ~/public/imasdb/${database}/4/111111/${run_number}                      || exit 1
      rm "IMAS_RUN_OUT" "IMAS_DB_OUT"
  fi
}

# --- Compare the results of the test case to the reference solution
function compare_results () {
    sed '/^\s*\*/d' jorek_namelist > file1
    sed '/^\s*\*/d' jorek_namelist_ref > file2  
  diff file1 file2                                             || exit 1
}
