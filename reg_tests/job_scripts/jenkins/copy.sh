#!/bin/bash
if [ $# -ne 2 ]; then
  echo "Usage: $0 <original job file> <new job file>"
  echo "       Copy the original file into the new job file, and "
  echo "       replace occcurences of original name into the copied file"
  exit 1
fi
ORIG=$1
TG=$2
PORIG="${ORIG%.*}"
PTG="${TG%.*}"
printf "Creating %-32s based on contents of the file $ORIG\n" $TG
sed "s/${PORIG}/${PTG}/g" $ORIG > $TG