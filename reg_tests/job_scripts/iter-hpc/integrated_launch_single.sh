#!/bin/bash
startdir=`readlink -f $(dirname $0)`
nonregdir=`readlink -f ${startdir}/../..` 
cd $nonregdir || exit 1
source job_scripts/iter-hpc/env.sh || exit 1
cd job_scripts/iter-hpc || exit 1
./$1.job || exit 1