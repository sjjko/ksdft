#!/bin/bash

#this script
#		runs make
#		makes soft links to scripts

#invoke with run.sh caseName

echo "run make " | tee log.setup

make release  | tee -a log.setup
ln ./bin/Release/ksdft++ release_ksdft++  | tee -a log.setup
ln -s script/cleanup.sh  | tee -a log.setup
ln -s script/run.sh  | tee -a log.setup
ln -s script/plotScript.sh  | tee -a log.setup
