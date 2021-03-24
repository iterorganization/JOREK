#!/bin/bash
LIST=`cd ../../..;./non_regression_tests/run_test.sh -L`
echo "== Creating job files =="
echo "   for Projet" $PROJ
for NAME in $LIST; do 
  TG=${NAME}.job; 
  sed -e "s/_PROJ_/${PROJ}/" -e "s/_CASE_/$NAME/g" case.skel > $TG; echo $TG;
  chmod +x $TG
done

