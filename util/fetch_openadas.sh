#!/usr/bin/env bash
# Retrieve the OPEN-ADAS ADF11 datafiles for a specific symbol.
# Downloads all files in parallel, removes the ones which were
# not found (the web service does not give a 404 but bad data)
#
# Author: Daan van Vugt <the JOREK team>
# Date: 2018-03-28

set -o pipefail
set -o errtrace
set -o nounset
set -o errexit

declare -a ADF11_sets=(acd scd qcd xcd ccd plt prb pls prc)

function usage() {
  echo "Retrieve OPEN-ADAS ADF11 datafiles from the OPEN-ADAS website"
  echo "and store them in files (type_\$symbol.dat) in the current directory."
  echo "by default we try to download all possible sets in ADF11 with this suffix,"
  echo "which are $(tr '[:lower:]' '[:upper:]' <<< ${ADF11_sets[@]})."
  echo "See http://open.adas.ac.uk/adf11 for more information."
  echo ""
  echo "Usage: $0 suffix [suffix2] [...]"
  echo "  where suffix is an identifier such as 50_w. See the OPEN-ADAS website"
  echo "  for more information on the different datasets."
  echo ""
  echo "Examples:"
  echo "  $0 50_w"
  echo ""
  echo "Options:"
  echo "  -h, --help: show this message"
}

if [ "$#" -le 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage && exit
fi

for element in "$@"; do
  type=(${element//_/ }) # an array of 50 and w
  for set in "${ADF11_sets[@]}"; do
    if [ ! -f $set$element.dat ]; then
      (wget -q http://open.adas.ac.uk/download/adf11/$set$type/$set$element.dat -O $set$element.dat &)
      echo "Downloading $set$element.dat"
    fi
  done
done
wait # for the downloader threads above
sleep 0.5 # wait for FS update
# Now check for bad files
for element in "$@"; do
  type=(${element//_/ }) # an array of 50 and w
  for set in "${ADF11_sets[@]}"; do
    if grep -q "You have an error" $set$element.dat; then
      rm $set$element.dat
      echo "Nonexisting $set$element.dat, skipping"
    fi
  done
done
