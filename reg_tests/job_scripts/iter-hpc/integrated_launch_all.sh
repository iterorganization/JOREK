#!/bin/bash
startdir=`readlink -f $(dirname $0)`
nonregdir=`readlink -f ${startdir}/../..` 
cd $nonregdir || exit 1
source job_scripts/iter-hpc/env.sh || exit 1
./launch_all.sh 