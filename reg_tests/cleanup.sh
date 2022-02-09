#!/bin/bash

startdir=`readlink -f $(dirname $0)`
codedir=`readlink -f ${startdir}/..` # Assumption about source code location
cd $codedir
for i in `ls -1d non_regression_tests/testcases/*`; do
  if [ -d $i ]; then
    cd $i
    echo "Cleaning up '$i'"
    
    source settings.sh
    
    rm -f *.h5 *.rst jorek_model* rst_*bin* *.tgz *.bck ${extra_remote_files}
    
    cd $codedir
  fi
done
