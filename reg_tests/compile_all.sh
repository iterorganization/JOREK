#!/bin/bash

startdir=`readlink -f $(dirname $0)`
codedir=`readlink -f ${startdir}/..` # Assumption about source code location
for dirname in ${startdir}/testcases/*; do
 if [ \( -d $dirname \)  -a \( ! -f $dirname/BROKEN \) ]; then
   name=$(basename $dirname)
   echo ""
   echo "===== Compile $name ====="
   (cd ${codedir}; non_regression_tests/run_test.sh $@ -p ${name} 1> log_$name.out 2> log_$name.err)  || (tail ${codedir}/log_$name.*; exit 1) || exit 1
 fi
done

echo ""
echo "COMPLETED: All test cases have been compiled successfully."

exit 0
