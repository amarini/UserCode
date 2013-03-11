#!/bin/bash
FILENAME=$1
ORIGIN="/eos/cms/store/user/webermat/zjgfiles_V00-07/"
DEST="/eos/cms/store/user/amarini/zjets_V00-07/"
EXEDIR=/afs/cern.ch/work/a/amarini/CMSSW_5_3_6/src/amarini/VPlusJets/post-process

mkdir -p /afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/

#---------- This script is executed on the remote system ------
cat >/afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.sh <<EOF
#!/bin/bash
export SCRAM_ARCH=slc5_amd64_gcc434 ; 
cd /afs/cern.ch/work/a/amarini/CMSSW_5_3_6/src ; eval \`scramv1 runtime -sh\` ; cd - ; 

cd \$WORKDIR
alias eos='/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select'
eos cp $ORGIN/$FILENAME ./
cd $EXECDIR
./PostProcess \$WORKDIR/$FILENAME 0
cd \$WORKDIR
eos cp $FILENAME $DEST

EOF

#---------- This script is executed on the local system ------
QUEUE=1nd
LOG=/afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.log

chmod a+rx /afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.sh 

echo " bsub -q ${QUEUE} -o ${LOG} -J ${FILENAME%%.root} < /afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.sh "
