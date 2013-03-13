#!/bin/bash
DEST="/eos/cms/store/user/amarini/zjets_V00-07"
EXEDIR=/afs/cern.ch/work/a/amarini/CMSSW_5_3_6/src/amarini/VPlusJetsAnalysis
SCRIPTDIR=/afs/cern.ch/work/a/amarini/VPlusJets_V00-07/Batch_scripts

mkdir -p $SCRIPTDIR

for i in {1..7} ;do
for config in data/config_MM.ini data/config_EE.ini data/config_TT.ini; do
#---------- This script is executed on the remote system ------
NAME=${config%%.ini}
NAME=${NAME##*/}_${i}
cat >${SCRIPTDIR}/script_$NAME.sh <<EOF
#!/bin/bash
export SCRAM_ARCH=slc5_amd64_gcc462;
cd /afs/cern.ch/work/a/amarini/CMSSW_5_3_6/src ; eval \`scramv1 runtime -sh\` ; cd - ; 

#mount eos
eosmount \${HOME}/eos
cd $EXEDIR

./ValidationPlot ${config} ${i}

#umount eos
eosumount \${HOME}/eos

EOF

#---------- This script is executed on the local system ------
QUEUE=8nh
LOG=${SCRIPTDIR}/log_$NAME.log

chmod a+rx ${SCRIPTDIR}/script_$NAME.sh 

[ -f "${LOG}" ] && rm ${LOG}
echo " bsub -q ${QUEUE} -o ${LOG} -J ${NAME} < ${SCRIPTDIR}/script_${NAME}.sh "

#-----------------------END------------------------------------
done
done
