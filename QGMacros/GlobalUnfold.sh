
#!/bin/bash
#root -l -b -q 'ComputeMixture.C("../Omog_DiJet_QCD_HT_Summer11.root","../Omog_QGStudies_*_Summer11.root",100,150,10,13,"Jet0","Jet0","omog")' 2>/dev/null | grep 'q/(q+g)'

#root 'Unfold.C("../Omog_DiJet_HT_Run2011_FULL.root","../Omog_QGStudies_Photon_Run2011_FULL.root","nCharged","(50,0,50)",100,150,10,13,.4,.85,false,"../Omog_DiJet_QCD_HT_Summer11.root","Jet0","Jet0","","omog")'

if [ "$1" == "" ]; then

export PTMIN=101
export PTMAX=127
export RHOMIN=10
export RHOMAX=13

else
export  PTMIN=$1
export  PTMAX=$2
export RHOMIN=$3
export RHOMAX=$4
fi;

if [ "$5" == "" ]; then
VARNAME="nCharged"
else
VARNAME=$5
fi;

if [ "$VARNAME" == "ptD" ];
then
RANGE="(50,0,1)"
fi;
if [ "$VARNAME" == "QGLikelihood" ]; then
RANGE="(50,0,1.0001)"
fi;
if [ "$VARNAME" == "nCharged" ]; then
RANGE="(50,0,50)"
fi;
if [ "$VARNAME" == "nNeutral" ]; then
RANGE="(50,0,50)"
fi;


##DCAP
#DIJETMC='dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/pandolf/Omog/Omog_DiJet_QCD_HT_Summer11.root'
#DIJETDATA='dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/pandolf/Omog/Omog_DiJet_HT_Run2011_FULL.root'
#PHJETMC='dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/pandolf/Omog/Omog_QGStudies_*_Summer11.root'
#PHJETDATA='dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/pandolf/Omog/Omog_QGStudies_Photon_Run2011_FULL.root'

##OMOG
DIR='/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini'
DIJETMC='/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_QCD_HT_Summer11.root'
DIJETDATA='/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_HT_Run2011_FULL.root'
PHJETMC='/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_*_Summer11.root'
PHJETDATA='/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_Photon_Run2011_FULL.root'

echo "Executing root - outputfile name=../Inversion/${VARNAME}_pt${PTMIN}_${PTMAX}_rho${RHOMIN}_${RHOMAX}.root"
echo "Range = ${RANGE}"

root  -l -b <<EOF
double q1,eq1,q2,eq2;
.L new/ComputeMixture.C+
.L Unfold.C
ComputeMixture("${DIJETMC}",${PTMIN},${PTMAX},${RHOMIN},${RHOMAX},"omog",&q1,&eq1);
ComputeMixture("${PHJETMC}",${PTMIN},${PTMAX},${RHOMIN},${RHOMAX},"omog",&q2,&eq2);
Unfold("${DIJETDATA}","${PHJETDATA}","${VARNAME}","${RANGE}",$PTMIN,$PTMAX,$RHOMIN,$RHOMAX,q1,q2,false,"${DIJETMC}","Jet0","Jet0","../Inversion/${VARNAME}_pt${PTMIN}_${PTMAX}_rho${RHOMIN}_${RHOMAX}.root","omog");

.x CompareSpectra.C(${PTMIN},${PTMAX},${RHOMIN},${RHOMAX},"../Inversion/PtJet0_pt${PTMIN}_${PTMAX}_rho${RHOMIN}_${RHOMAX}.pdf");
.q
EOF

echo -e '\033[01;32m DONE \033[00m'
