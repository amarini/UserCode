//Uptadet To omog:
int PlotChiSquare2(const char*varName="nChargedJet0",char what='D',int PtMin=160,int PtMax=200,int RhoMin=4,int RhoMax=6,const char *fileName=""){
TChain *a=new TChain("omog");
TChain *b=new TChain("omog");
TChain *a1=new TChain("omog");
TChain *b1=new TChain("omog");
if(what=='D'){
printf("\033\[01;31mDijet\033\[0m\n");
a->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_HT_Run2011_FULL.root");
b->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_QCD_HT_Summer11.root");
a1->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_DiJet_HT_Run2011_FULL.root");
b1->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_DiJet_QCD_HT_Summer11.root");
}else if(what=='P'){
printf("\033\[01;31mPhotons\033\[0m\n");
a->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_Photon_Run2011_FULL.root");
b->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_*_Summer11.root");
a1->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_QGStudies_Photon_Run2011_FULL.root");
b1->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_QGStudies_*_Summer11.root");
}
a->AddFriend(a1);
b->AddFriend(b1);
TH1F*Q=new TH1F("nC_q","nC_q",100,-0.5,99.5);
TH1F*G=new TH1F("nC_g","nC_g",100,-0.5,99.5);
TH1F*T=new TH1F("nC_t","nC_t",100,-0.5,99.5);
TH1F*T2=new TH1F("nC_t2","nC_t2",100,-0.5,99.5);
Q->Sumw2();G->Sumw2();T->Sumw2();T2->Sumw2();
char name[1023];
char sel[1023];
sprintf(name,"%s>>nC_q",varName);
sprintf(sel,"eventWeight*(ptJet0>%d && ptJet0<%d && %d<rhoPF && rhoPF<%d && abs(pdgIdJet0)<4 && pdgIdJet0!=0)",PtMin,PtMax,RhoMin,RhoMax);
b->Draw(name,sel);
sprintf(name,"%s>>nC_g",varName);
sprintf(sel,"eventWeight*(ptJet0>%d && ptJet0<%d && %d<rhoPF && rhoPF<%d && pdgIdJet0==21)",PtMin,PtMax,RhoMin,RhoMax) ;
b->Draw(name,sel);
sprintf(name,"%s>>nC_t2",varName);
sprintf(sel,"eventWeight*(ptJet0>%d && ptJet0<%d && %d<rhoPF && rhoPF<%d )",PtMin,PtMax,RhoMin,RhoMax);
b->Draw(name,sel);
sprintf(name,"%s>>nC_t",varName);
sprintf(sel,"(ptJet0>%d && ptJet0<%d && %d<rhoPF && rhoPF<%d )",PtMin,PtMax,RhoMin,RhoMax);
a->Draw(name,sel);

//.L Macros/new/MinChiSquare.C 
gSystem->Load("new/MinChiSquare_C.so");
TGraph *graph=new TGraph();graph->SetName("data");
MinChiSquare::GetMinChiSquare(Q,G,T,0.01,"WW CHI2",graph);
TGraph *graph2=new TGraph();graph2->SetName("mc");
MinChiSquare::GetMinChiSquare(Q,G,T2,0.01,"WW CHI2",graph2);

graph->SetMarkerStyle(20);
graph->SetMarkerColor(kBlack);
graph2->SetMarkerStyle(24);
graph2->SetMarkerColor(kRed);

TCanvas *c1=new TCanvas("c1","c1");
graph->Draw("ALP");
graph2->Draw("P SAME");

if(fileName[0]!='\0')c1->SaveAs(fileName);
return 0;
}
