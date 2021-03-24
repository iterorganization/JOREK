#!/bin/bash

if [ -z "${DAV_URL}" ]; then
    DAV_URL="http://jorek.eu/dav_nrt"
fi

while [ $# -gt 0 ]; do 
  wget -q --user=nrt --password=nrt_21745XtL ${DAV_URL}/$1
  if [ $? -eq 0 ]; then
    echo "successfully downloaded '$1'"
    shift
  else
    echo "error downloading '$1'"
    exit 1
  fi
done

