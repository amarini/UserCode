
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
else 
RANGE="(50,0,50)"
fi;

##VAR=$( eval root -l -b -q \'ComputeMixture.C\(\"../Omog_DiJet_QCD_HT_Summer11.root\",\"../Omog_QGStudies_\*_Summer11.root\",$PTMIN,$PTMAX,$RHOMIN,$RHOMAX,\"Jet0\",\"Jet0\",\"omog\"\)\' 2>&1 | grep 'q/q+g' | cut -d ' ' -f 3 | cut -d'=' -f 2 | tr '\n' ' ')
##DI=$(echo $VAR| cut -d ' ' -f1)
##PH=$(echo $VAR| cut -d ' ' -f2)
##
##echo $DI $PH
## eval root -l -q -b \'Unfold.C\(\"../Omog_DiJet_HT_Run2011_FULL.root\",\"../Omog_QGStudies_Photon_Run2011_FULL.root\",\"$VARNAME\",\"${RANGE}\",$PTMIN,$PTMAX,$RHOMIN,$RHOMAX,$DI,$PH,false,\"../Omog_DiJet_QCD_HT_Summer11.root\",\"Jet0\",\"Jet0\",\"../Inversion/${VARNAME}_pt${PTMIN}_${PTMAX}_rho${RHOMIN}_${RHOMAX}.root\",\"omog\"\)\'  
##

DIJETMC='../Omog_DiJet_QCD_HT_Summer11.root'
DIJETDATA='../Omog_DiJet_HT_Run2011_FULL.root'
PHJETMC='../Omog_QGStudies_*_Summer11.root'
PHJETDATA='../Omog_QGStudies_Photon_Run2011_FULL.root'
echo "Executing root - outputfile name=../Inversion/${VARNAME}_pt${PTMIN}_${PTMAX}_rho${RHOMIN}_${RHOMAX}.root"


root  -l -b <<EOF
double q1,eq1,q2,eq2;
.L new/ComputeMixture.C+
.L Unfold.C
ComputeMixture("${DIJETMC}",${PTMIN},${PTMAX},${RHOMIN},${RHOMAX},"omog",&q1,&eq1);
ComputeMixture("${PHJETMC}",${PTMIN},${PTMAX},${RHOMIN},${RHOMAX},"omog",&q2,&eq2);
Unfold("${DIJETDATA}","${PHJETDATA}","${VARNAME}","${RANGE}",$PTMIN,$PTMAX,$RHOMIN,$RHOMAX,q1,q2,false,"${DIJETMC}","Jet0","Jet0","../Inversion/${VARNAME}_pt${PTMIN}_${PTMAX}_rho${RHOMIN}_${RHOMAX}.root","omog");
.q
EOF

