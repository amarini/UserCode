#!/bin/bash

for config in data/config_MM.ini data/config_EE.ini ; do
for type in {0..1} ;do
#
#for j in {0..9} ; do

	
bsub -q 1nh -o output/log_QG.txt <<EOF
export "SCRAM_ARCH=slc5_amd64_gcc462"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
cd ~/work/CMSSW_5_3_6/src/amarini/VPlusJetsAnalysis/CommonFramework
cmsenv
./qgexe data/config_LLPt.ini $config $type
EOF

done
done
#done

