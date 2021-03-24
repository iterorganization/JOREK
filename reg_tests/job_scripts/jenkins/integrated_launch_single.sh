#!/bin/bash
startdir=`readlink -f $(dirname $0)`
nonregdir=`readlink -f ${startdir}/../..` 
cd $nonregdir || exit 1
source job_scripts/jenkins/env.sh || exit 1
cd job_scripts/jenkins || exit 1
./$1.job || exit 1