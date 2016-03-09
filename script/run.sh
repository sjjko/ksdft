#!/bin/bash

#this script
#		cleans up the postprocessing files
#		runs ksdft++ program
#		does postprocessing
#		moves images and data to POSTPROCESSING/EPS and POSTPROCESSING/DATA folders
#		generates a 

#invoke with run.sh caseName

if [ -z "$1" ]; then
  echo "supply caseName as first argument to this script" | tee log.runsh
  exit
fi

CASENAME=$1
echo "run ksdft with case " $CASENAME | tee -a log.runsh
echo "====================" | tee -a log.runsh
echo "call cleanup script" | tee -a log.runsh
bash cleanup.sh | tee -a log.runsh
echo "finished cleanup" | tee -a log.runsh
echo "====================" | tee -a log.runsh
echo "call ksdft++" | tee -a log.runsh
rm ./log.minidft
./release_ksdft++ $CASENAME | tee -a log.minidft log.runsh
echo "finished ksdft++" | tee -a log.runsh
echo "====================" | tee -a log.runsh
echo "call post processing routine" | tee -a log.runsh
bash plotScript.sh $CASENAME  | tee -a log.run | tee -a log.runsh
echo "finished postprocessing routine" | tee -a log.runsh
echo "cp latex pdf to ./doc folder" | tee -a log.runsh
cp latex.pdf ./doc/ksdft++.pdf
echo "====================" | tee -a log.runsh
