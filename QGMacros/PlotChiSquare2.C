
int PlotChiSquare2(const char*varName="nChargedJet0"){
TChain *a=new TChain("omog");
TChain *b=new TChain("omog");
a->Add("Omog_DiJet_HT_Run2011_FULL.root");
b->Add("Omog_DiJet_QCD_HT_Summer11.root");
TH1F*Q=new TH1F("nC_q","nC_q",100,-0.5,99.5);
TH1F*G=new TH1F("nC_g","nC_g",100,-0.5,99.5);
TH1F*T=new TH1F("nC_t","nC_t",100,-0.5,99.5);
TH1F*T2=new TH1F("nC_t2","nC_t2",100,-0.5,99.5);
Q->Sumw2();G->Sumw2();T->Sumw2();T2->Sumw2();
char name[1023];
sprintf(name,"%s>>nC_q",varName);
b->Draw(name,"eventWeight*(ptJet0>160 && ptJet0<200 && 4<rhoPF && rhoPF<6 && abs(pdgIdJet0)<4 && pdgIdJet0!=0)");
sprintf(name,"%s>>nC_g",varName);
b->Draw(name,"eventWeight*(ptJet0>160 && ptJet0<200 && 4<rhoPF && rhoPF<6 && pdgIdJet0==21)") ;
sprintf(name,"%s>>nC_t2",varName);
b->Draw(name,"eventWeight*(ptJet0>160 && ptJet0<200 && 4<rhoPF && rhoPF<6 )");
sprintf(name,"%s>>nC_t",varName);
a->Draw(name,"(ptJet0>160 && ptJet0<200 && 4<rhoPF && rhoPF<6 )");

//.L Macros/new/MinChiSquare.C 
gSystem->Load("Macros/new/MinChiSquare_C.so");
TGraph *graph=new TGraph();graph->SetName("data");
MinChiSquare::GetMinChiSquare(Q,G,T,0.01,"WW CHI2",graph);
TGraph *graph2=new TGraph();graph2->SetName("mc");
MinChiSquare::GetMinChiSquare(Q,G,T2,0.01,"WW CHI2",graph2);

graph->SetMarkerStyle(20);
graph->SetMarkerColor(kBlack);
graph2->SetMarkerStyle(24);
graph2->SetMarkerColor(kRed);

graph->Draw("ALP");
graph2->Draw("P SAME");

return 0;
}
