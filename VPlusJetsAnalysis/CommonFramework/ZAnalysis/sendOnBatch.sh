#!/bin/bash

#./zexe data/config_MM.ini 0 0 100

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
for j in {0..9} ; do

if [ "$1" == "check" ] ; then

[ "$config" == "data/config_MM.ini" ] && CHID=4;
[ "$config" == "data/config_EE.ini" ] && CHID=1;
[ "$config" == "data/config_TT.ini" ] && CHID=2;

RESUBMIT=0

[ -f $OUTDIR/output_type${type}_chid${CHID}_${j}_10.root ] || { echo "NotFinished:" ./zexe $config $type $j 10 ; RESUBMIT=1; }
[ -f $OUTDIR/output_type${type}_chid${CHID}_${j}_10.root ] && [ `ls -l ${OUTDIR}/output_type${type}_chid${CHID}_${j}_10.root | tr -s ' ' | cut -d' ' -f5` -le 100 ] && echo "NotFinished:" ./zexe $config $type $j 10  && RESUBMIT=1;

[ "$RESUBMIT" == "1" ] && [ "$2" == "resubmit" ] && bsub -q 8nh -o output/log.txt <<ABC
export "SCRAM_ARCH=slc5_amd64_gcc462"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
cd ~/work/CMSSW_5_3_6/src/amarini/VPlusJetsAnalysis/CommonFramework
cmsenv
./zexe $config $type $j 10
ABC

else
if [ "$1" != "add" ]; then
	
bsub -q 1nh -o output/log.txt <<EOF
export "SCRAM_ARCH=slc5_amd64_gcc462"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
cd ~/work/CMSSW_5_3_6/src/amarini/VPlusJetsAnalysis/CommonFramework
cmsenv
./zexe $config $type $j 10
EOF
fi;
fi;
done
done
done

if [ "$1" == "add" ]; then
for CHID in 1 2 4 ; do
for type in {0..6} ;do
hadd $OUTDIR/output_type${type}_chid${CHID}.root $OUTDIR/output_type${type}_chid${CHID}_*_10.root 
done; done;
fi;
