#!/bin/bash

startdir=`readlink -f $(dirname $0)`
codedir=`readlink -f ${startdir}/..` # Assumption about source code location
for dirname in ${startdir}/testcases/*; do
  if [ -d $dirname ] && [ ! -e "$dirname/BROKEN" ]; then
    name=$(basename $dirname)
    echo ""
    echo "===== Run $name ====="
    
    cd ${codedir}
    non_regression_tests/run_test.sh $@ -n ${name} 1> log_$name.out 2> log_$name.err
    ierr=$?
    
    cat ${codedir}/log_$name.out | sed -e '/^ *$/d' | tail -n 1
    if [ ! $ierr -eq 0 ]; then
      cat ${codedir}/log_$name.err | sed -e '/^ *$/d' | tail -n 1
      echo ""
      echo "ERROR running testcase $name! Details in the logfiles:"
      echo "  ${codedir}/log_$name.out"
      echo "  ${codedir}/log_$name.err"
      exit 1
    fi
  fi
done
echo ""
echo "COMPLETED: All test cases have been carried out successfully."
exit 0
