#!/bin/bash

# This script upload a jorek test case.
TESTDIR="$1"
if [ -z "${DAV_URL}" ]; then
    DAV_URL="http://jorek.eu/dav_nrt/"
fi
if [ ! -f ${TESTDIR}/settings.sh ]; then
  printf "This test name does not exist\n"
  exit 1
fi

cd  ${TESTDIR}
testcasedir=`readlink -f ${PWD}`
source ./settings.sh

VERSION="_`cat end.h5 $extra_remote_files | md5sum | sed -e 's/^ //' -e 's/ .*$//'`"
echo ${VERSION} > .version

printf "date: "      > .info
date                >> .info
printf "hostname: " >> .info
hostname -f         >> .info
printf "branch: "   >> .info
git rev-parse --abbrev-ref HEAD >> .info
printf "user: "     >> .info
echo $USER          >> .info

testname=$(basename $testcasedir)
TGZFILE=${testname}${VERSION}.tgz

echo "Creating tarball ${TGZFILE}"
tar cvzf ${TGZFILE} .info begin.h5 end.h5 ${extra_remote_files} || exit 1
 
echo "Uploading ${TGZFILE}"
TESTNAME=$(basename $TESTDIR)
curl -s -u nrt:nrt_21745XtL -T ${TGZFILE} ${DAV_URL}
if [ $? -eq 0 ]; then
  printf "Success\n"
  exit 0
else
  printf "Failed\n"
  exit 1
fi
