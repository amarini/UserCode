#!/bin/bash

#./zexe data/config_MM.ini 0 0 100


for config in data/config_MM.ini data/config_EE.ini data/config_TT.ini; do
for type in {0..6} ;do
#
for j in {0..9} ; do

if [ "$1" == "check" ] ; then

[ "$config" == "data/config_MM.ini" ] && CHID=4;
[ "$config" == "data/config_EE.ini" ] && CHID=1;
[ "$config" == "data/config_TT.ini" ] && CHID=2;

[ -f output/output_type${type}_chid${CHID}_${j}_10.root ] || echo "NotFinished:" ./zexe $config $type $j 10 ;
[ -f output/output_type${type}_chid${CHID}_${j}_10.root ] && [ `ls -l output/output_type${type}_chid${CHID}_${j}_10.root | tr -s ' ' | cut -d' ' -f5` -le 100 ] && echo "NotFinished:" ./zexe $config $type $j 10 ; 

else

bsub -q 1nh -o output/log.txt <<EOF
export "SCRAM_ARCH=slc5_amd64_gcc462"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
cd ~/work/CMSSW_5_3_6/src/amarini/VPlusJetsAnalysis/CommonFramework
cmsenv
./zexe $config $type $j 10
EOF

fi;
done
done
done
