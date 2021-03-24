#!/bin/bash
LIST=`cd ../../..;./non_regression_tests/run_test.sh -L`
echo "== Creating job files =="
for NAME in $LIST; do 
  TG=${NAME}.job; 
  sed  "s/_CASE_/$NAME/g" case.skel > $TG; echo $TG; 
  chmod +x $TG
done

