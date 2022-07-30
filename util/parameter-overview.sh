#!/bin/bash

# 
# This script extracts lists of input parameters from all models and creates an overview that can be copied to our wiki.
# This way, comments on the input parameters in the code are directly reflected in the wiki.
# (c) Matthias Hoelzl, 2019
# 

startdir=`pwd`

outfile="parameter-overview.txt"

models=""
for i in `ls -1d models/model*/initialise_parameters.f90`; do
  model=`echo $i | sed -e 's|models/model||' -e 's|/.*||'`
  if [ $model = "001" ] || [ $model = "002" ] || [ $model = "003" ] || [ $model = "004" ] || [ $model = "005" ] || [ $model = "305" ] || [ $model = "306" ] || [ $model = "333" ]; then
    continue
  fi
  models="$models $model"
  grep /in1/ models/model$model/initialise_parameters.f90 -A 999 | grep -v "^#" | grep "&" -A 1 | tr -d '\n' | sed -e 's/[&,]/ /g' -e 's/\t/ /g' -e 's/  */ /g' -e 's/,//g' -e 's|^.*/ ||' | tr ' ' '\n' | sort | uniq > tmp_${model}_$$
done

cat tmp_*_$$ | sort | uniq > tmp_$$

function insert_header() {
  echo -n "^  parameter  ^  default  ^  description  ^  " >> $outfile
  for model in $models; do
    echo -n "$model  ^  " >> $outfile
  done
  echo "" >> $outfile
}

function find_variables() {
  file=$1
  grep -B 9999 'contains' $file | egrep '^[^!]*::' | sed -e 's/^.*:://' -e s'/!.*$//' -e 's/([^)]*)//g' -e 's/=.*//' -e 's/&//' -e 's/,/ /g' -e 's/$/ /' | tr -d '\n' | sed -e 's/\t/ /g' -e 's/  */ /g' | tr ' ' '\n'
}

function is_in_namelists() {
  param=$1
  anymatch=0
  for model in $models; do
    matches=`egrep -i -x $param tmp_${model}_$$ | wc -l`
    if [ "$matches" == "1" ]; then
      anymatch=1
    fi
  done
  echo $anymatch
}

cat communication/broadcast_phys.f90 > tmp_$$_comm
egrep -A 9999 "^ *subroutine broadcast_vacuum" vacuum/vacuum.f90 | grep -B 9999 "end subroutine broadcast_vacuum" >> tmp_$$_comm

function is_communicated() {
  param=$1
  num_found=`egrep -i "^[^!]*\([ \t]*$param[ ,(]"  tmp_$$_comm | wc -l`
  if [ "$num_found" == "0" ]; then
    echo "Warning: Parameter $param might not get communicated!"
  fi
}

find_variables models/phys_module.f90 > tmp_$$_phys
find_variables vacuum/vacuum.f90      > tmp_$$_vacu

rm -f $outfile

echo "===== phys_module =====" >> $outfile
insert_header
k=0
for param in `cat tmp_$$_phys`; do
  if [ "`is_in_namelists $param`" == "1" ]; then
    k=$[k+1]
    m=$(( $k % 20 ))
    if [ $m -eq 0 ]; then
      insert_header
    fi
    is_communicated $param
    description=`egrep -i "^[^!]* $param[( ]" models/phys_module.f90 | grep "!" | sed -e 's/^.*![< ]*//' -e 's/\\\f//g' | tr '\n' ';' | sed -e 's/;$//' -e 's/((/( (/' -e 's/))/) )/' -e s'|//|/ /|'`
    default=`egrep -i "^[^!] $param[ =(]" models/preset_parameters.f90 | sed -e 's/^[^!]*= *//' -e 's/!.*$//' -e 's| *(/ *||' -e 's| */) *||' -e 's/d0//g' -e 's/rst_hdf5_version_supported//' -e 's/[ \t]*$//' | tr '\n' ' '`
    echo -n "| **$param** | $default | $description |" >> $outfile
    for model in $models; do
      matches=`egrep -x $param tmp_${model}_$$ | wc -l | sed -e 's/0/ /' -e "s/1/x/"`
      echo -n "  $matches  |" >> $outfile
    done
    echo "" >> $outfile
  fi
done

echo "===== vacuum =====" >> $outfile
insert_header
k=0
for param in `cat tmp_$$_vacu`; do
  if [ "`is_in_namelists $param`" == "1" ]; then
    k=$[k+1]
    m=$(( $k % 20 ))
    if [ $m -eq 0 ]; then
      insert_header
    fi
    is_communicated $param
    description=`egrep -i "^[^!]* $param[( ]" vacuum/vacuum.f90 | grep "!" | sed -e 's/^.*![< ]*//' -e 's/\\\f//g' | tr '\n' ';' | sed -e 's/;$//' -e 's/((/( (/' -e 's/))/) )/' -e s'|//|/ /|'`
    default=`egrep -A 9999 "^ *subroutine vacuum_preset" vacuum/vacuum.f90 | grep -B 9999 "end subroutine vacuum_preset" | grep "$param[ =(]" | sed -e 's/^[^!]*= *//' -e 's/!.*$//' -e 's/d0//g'`
    echo -n "| **$param** | $default | $description |" >> $outfile
    for model in $models; do
      matches=`egrep -x $param tmp_${model}_$$ | wc -l | sed -e 's/0/ /' -e "s/1/x/"`
      echo -n "  $matches  |" >> $outfile
    done
    echo "" >> $outfile
  fi
done

rm -f tmp*_$$*
