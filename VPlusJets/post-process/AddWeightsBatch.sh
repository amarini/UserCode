#!/bin/bash
FILENAME=$1
ORIGIN="/eos/cms/store/user/webermat/zjgfiles_V00-07"
DEST="/eos/cms/store/user/amarini/zjets_V00-07"
EXEDIR=/afs/cern.ch/work/a/amarini/CMSSW_5_3_6/src/amarini/VPlusJets/post-process
EOS='/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select'

mkdir -p /afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/

echo $FILENAME | grep "PartI\.root" >/dev/null 2>/dev/null
MATCHPARTI=$?
echo $FILENAME | grep "Part.*\.root" >/dev/null 2>/dev/null
MATCHPARTX=$?

#---------- This script is executed on the remote system ------
cat >/afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.sh <<EOF
#!/bin/bash
export SCRAM_ARCH=slc5_amd64_gcc462;
cd /afs/cern.ch/work/a/amarini/CMSSW_5_3_6/src ; eval \`scramv1 runtime -sh\` ; cd - ; 

cd \$WORKDIR


FILE=$FILENAME

[ ${MATCHPARTI} == 0 ] && { $EOS cp -a $ORIGIN/${FILENAME%%Part*}Part* ./ ; FILE=${FILENAME%%_Part*}.root; hadd \$FILE ${FILENAME%%Part*}Part*.root; }
[ ${MATCHPARTI} != 0 ] &&  [ ${MATCHPARTX} == 0 ] && exit 0;

[ ${MATCHPARTI} != 0 ] && { ${EOS} cp -a ${ORIGIN}/${FILENAME} ./ ;}

cd $EXECDIR
./PostProcess \$WORKDIR/\$FILE 0
cd \$WORKDIR
${EOS} cp -a $FILE ${DEST}/

EOF

#---------- This script is executed on the local system ------
QUEUE=1nd
LOG=/afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.log

chmod a+rx /afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.sh 

echo " bsub -q ${QUEUE} -o ${LOG} -J ${FILENAME%%.root} < /afs/cern.ch/work/a/amarini/VPlusJets_PostProcess/${FILENAME%%.root}.sh "

#alias eos='/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select'
