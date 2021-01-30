#!/bin/bash

# 
# This script checks whether all input parameters are properly communicated.
# (c) Matthias Hoelzl, 2021
# 

startdir=`pwd`

if [ $# -ne 0 ]; then
  if [ $# -eq 1 ]; then
    if [ $1 == "-h" ]; then
      echo "Usage: ./util/check-communicated.sh [-h] [-o]"
      echo "-h print help"
      echo "-o only check whether parameters communicated"
      exit
    elif [ $1 == "-o" ]; then
      onlycheck="1"
    else
      echo "Wrong parameters"
      exit 1
    fi
  else
    echo "Wrong number of parameters"
    exit 1
  fi
fi

models=""
for i in `ls -1d models/model*/initialise_parameters.f90`; do
  model=`echo $i | sed -e 's|models/model||' -e 's|/.*||'`
  if [ $model = "001" ] || [ $model = "305" ] || [ $model = "306" ] || [ $model = "555" ] || [ $model = "701" ]; then
    continue
  fi
  #echo "Checking model $model"
  models="$models $model"
  grep /in1/ models/model$model/initialise_parameters.f90 -A 999 | grep "&" -A 1 | tr -d '\n' | sed -e 's/[&,]/ /g' -e 's/\t/ /g' -e 's/  */ /g' -e 's/,//g' -e 's|^.*/ ||' | tr ' ' '\n' | sort | uniq > tmp_${model}_$$
  #echo "Parameter list for model $model created" 
done

cat tmp_*_$$ | sort | uniq > tmp_$$
#echo "Global parameter list created" 

cat communication/broadcast_phys.f90 | grep -i "call mpi" > tmp_$$_comm
egrep -A 9999 "^ *subroutine broadcast_vacuum" vacuum/vacuum.f90 | grep -B 9999 "end subroutine broadcast_vacuum" | grep -i "call mpi" >> tmp_$$_comm
communicationlist=`cat tmp_$$_comm`
#echo "Communication list created"

function is_communicated() {
  param=$1
  num_found=`echo "$communicationlist" | egrep -i "^[^!]*\([ \t]*$param[ ,(]" | wc -l`
  if [ "$num_found" == "0" ]; then
    echo "Warning: Parameter $param might not get communicated!"
    err=1
  fi
}

err=0
for param in `cat tmp_$$`; do
  echo "Checking parameter $param"
  is_communicated $param
done

rm -f tmp*_$$*

if [ $err -ne 0 ]; then
  echo ""
  echo "WARNING: Some parameters seem not to be communicated!"
  echo ""
fi

exit $err
