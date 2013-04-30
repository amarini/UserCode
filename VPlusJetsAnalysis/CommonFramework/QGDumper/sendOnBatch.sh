#!/bin/bash

#./qgexe data/config_MM.ini 0 0 100

echo --$1--
echo "--->$2<---"

[ "$1" == "check" ] && echo "Going to check ..."
[ "$1" == "check" ] && [ "$2" == "resubmit" ] && echo "  -->  And resubmit ..."
[ "$1" == "add" ] && echo "Going to add ..."

[ "$1" != "check" ] && [ "$1" != "add" ] && echo "Going to submit ..."

read

OUTDIR=~/work/VPlusJets_V00-07/Output_V00-07_CommonFramework
for config in data/config_MM.ini data/config_EE.ini data/config_TT.ini; do
for type in {0..6} ;do
#
#for j in {0..9} ; do

if [ "$1" == "check" ] ; then

[ "$config" == "data/config_MM.ini" ] && CHID=4;
[ "$config" == "data/config_EE.ini" ] && CHID=1;
[ "$config" == "data/config_TT.ini" ] && CHID=2;

RESUBMIT=0

	
bsub -q 1nh -o output/log_QG.txt <<EOF
export "SCRAM_ARCH=slc5_amd64_gcc462"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
cd ~/work/CMSSW_5_3_6/src/amarini/VPlusJetsAnalysis/CommonFramework
cmsenv
./qgexe $config $type
EOF
fi;
fi;
done
done
#done

