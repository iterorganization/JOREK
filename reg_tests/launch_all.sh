#!/bin/bash

startdir=`readlink -f $(dirname $0)`
codedir=`readlink -f ${startdir}/..` # Assumption about source code location
if [ -z "$JOREK_HOST" ]; then
    echo "JOREK_HOST environment variable is not defined, can not launch test cases"
    exit 1
fi
if [ -z "$BATCHCOMMAND" ]; then
    echo "BATCHCOMMAND environment variable is not defined, can not launch test cases"
    exit 1
fi
cd  || exit 1
for dirname in ${startdir}/testcases/*; do
 name=$(basename $dirname)
 if [ \( -d $dirname \) -a \( -f ${startdir}/job_scripts/${JOREK_HOST}/${name}.job \)  -a \( ! -f $dirname/BROKEN \) ]; then
   echo "== Launch job $name from directory job_scripts/${JOREK_HOST}"
   (cd ${startdir}/job_scripts/${JOREK_HOST}; ${BATCHCOMMAND} ${name}.job || exit 1)
 fi
done
exit 0
