#!/bin/bash

startdir=`readlink -f $(dirname $0)`
nonregdir=`readlink -f ${startdir}/../..` 
cd $nonregdir || exit 1
rm -f testcases/${1}/*tgz 2>/dev/null

set -v 
cp ../Make.inc/Makefile.iter_hpc_h5_linux_intel.inc ../Makefile.inc || exit 1
source job_scripts/iter-hpc/env.sh || exit 1
./get_all_data.sh || exit 1
cd ..
make cleanall
non_regression_tests/run_test.sh $1 -p

