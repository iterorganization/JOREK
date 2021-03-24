#!/bin/bash
source env.sh
LIST=`cd ../../..;./non_regression_tests/run_test.sh -L`
echo "== Creating job files =="
for NAME in $LIST; do 
  TG=${NAME}.job; 
  if [ "$TG" != "tear_circ_303.job" ]; then
    ./copy.sh tear_circ_303.job $TG;    
  fi
  chmod +x $TG
done
