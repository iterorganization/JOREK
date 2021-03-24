#!/bin/bash 

#output information about the regression before the test


#moved from env.sh 
module avail MUMPS
module avail SCOTCH
module avail PaStiX
module avail ParMETIS
module show PaStiX/5.2.3-intel-2018a
module show SCOTCH/6.0.4-intel-2018a
module show MUMPS/5.1.2-intel-2018a-metis
module show ParMETIS/4.0.3-intel-2018a
module list
ls $EBROOTPASTIX/lib
ls $EBROOTSCOTCH/lib
ls $EBROOTMUMPS/lib
ls $EBROOTPARMETIS/lib


echo -e "\n\n\n\n"
env
pwd
df -h
ls -lrt
ls -lrt non_regression_tests/testcases/
cat /proc/cpuinfo
git status

du -sh .
