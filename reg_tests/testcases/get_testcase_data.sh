#!/bin/bash

# This script download a jorek test case.
TESTNAME=`echo "$1" | sed -e 's|/$||'`
if [ -z "${DAV_URL}" ]; then
    DAV_URL="http://jorek.eu/dav_nrt"
fi
if [ ! -f ${TESTNAME}/settings.sh ]; then
  printf "This test name does not exist\n"
  exit 1
fi

cd  ${TESTNAME}

VERSION=""
if [ -f .version ]; then
  VERSION="`cat .version`"
fi

TGZFILE="${TESTNAME}${VERSION}.tgz"

if [ ! -f ${TGZFILE} ]; then 
  echo "Downloading ${TGZFILE}"
  wget -q --user=nrt --password=nrt_21745XtL ${DAV_URL}/${TGZFILE}
  returncode=$?
  if [ $returncode -ne 0 ]; then
    cat <<EOF
####################################################################
  Failed to automatically download from web site.
  Please download reference data ${TGZFILE} yourself 
  from http:/jorek.eu/dav_nrt and
  copy it into testcases/${TESTNAME} directory.
  Launch this script again to decompress the archive.
####################################################################
EOF
    exit 1
  fi
else 
    cat <<EOF
####################################################################
 No downloading performed because ${TGZFILE} already
 exists, remove it if you want to download:
   rm testcases/${TESTNAME}/${TGZFILE}
####################################################################
EOF
    exit 1
fi

echo "Uncompress ${TGZFILE}"
tar xvzf ${TGZFILE}
