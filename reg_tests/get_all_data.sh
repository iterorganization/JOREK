#!/bin/bash

startdir=$(dirname "$(readlink -f $0)")
codedir=`readlink -f "${startdir}/.."` # Assumption about source code location
cd "${startdir}/testcases" || exit 1
for name in *; do
  if [ -d $name ]; then
    echo ""
    echo "===== Getting data for $name ====="
    ./get_testcase_data.sh $name
  fi
done
exit 0
