//it takes Di or Ph+jets and Plot the QGFit4 with components
#include "/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/QGMacros/new/CMSLatex.C"
int MergePlot(char what='P',int ptMin=101,int ptMax=127,int rhoMin=4,int rhoMax=6,long long entries=1000000000)
{
char name[1023],sel[1023];

TChain *Di=new TChain("omog");
TChain *Di_QG=new TChain("omog");

TChain *Di_MC=new TChain("omog");
TChain *Di_QGMC=new TChain("omog");

//DATA
if(what=='P'){
Di->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_Photon_Run2011_FULL.root");
Di_QG->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_QGStudies_Photon_Run2011_FULL.root");
Di->AddFriend(Di_QG);
}else{
Di->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_HT_Run2011_FULL.root");
Di_QG->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_DiJet_HT_Run2011_FULL.root");
Di->AddFriend(Di_QG);
}
//MC - Use This for Ph
if(what=='P'){
Di_MC->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_*_Summer11.root");
Di_QGMC->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_QGStudies_*_Summer11.root");
Di_MC->AddFriend(Di_QGMC);
}else{
Di_MC->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_QCD_HT_Summer11.root");
Di_QGMC->Add("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Friend_DiJet_QCD_HT_Summer11.root");
Di_MC->AddFriend(Di_QGMC);
}

TH1F* QG4_Di=new TH1F("QG4_Di","QG4_Di",20,0,1);//data distributions

sprintf(sel,"eventWeight*(%d<ptJet0 &&ptJet0<%d &&%d<rhoPF &&rhoPF<%d && abs(etaJet0<2))",ptMin,ptMax,rhoMin,rhoMax);
Di->Draw("QGFit4>>QG4_Di",sel,"goff",entries);
QG4_Di->SetMarkerColor(kBlack);
QG4_Di->SetMarkerStyle(20);
QG4_Di->Scale(1./QG4_Di->Integral("width"));

TH1F* QG4_Q=new TH1F("QG4_Q","QG4_Q",20,0,1);//data distributions
TH1F* QG4_B=new TH1F("QG4_B","QG4_B",20,0,1);//data distributions
TH1F* QG4_G=new TH1F("QG4_G","QG4_G",20,0,1);
TH1F* QG4_U=new TH1F("QG4_U","QG4_U",20,0,1);//data distributions
QG4_Q->SetFillColor(kRed-4);    QG4_Q->SetLineColor(kBlack);QG4_Q->SetLineWidth(2); 
QG4_G->SetFillColor(kBlue-4);   QG4_G->SetLineColor(kBlack);QG4_G->SetLineWidth(2); 
QG4_B->SetFillColor(kYellow);   QG4_B->SetLineColor(kBlack);QG4_B->SetLineWidth(2); 
QG4_U->SetFillColor(kGray);     QG4_U->SetLineColor(kBlack);QG4_U->SetLineWidth(2); 
sprintf(sel,"eventWeight*(pdgIdJet0==21 &&%d<ptJet0 &&ptJet0<%d &&%d<rhoPF &&rhoPF<%d && abs(etaJet0<2))",ptMin,ptMax,rhoMin,rhoMax);
Di_MC->Draw("QGFit4>>QG4_G",sel,"goff",entries);
sprintf(sel,"eventWeight*(pdgIdJet0!=0 && abs(pdgIdJet0)<5 &&%d<ptJet0 &&ptJet0<%d &&%d<rhoPF &&rhoPF<%d && abs(etaJet0<2))",ptMin,ptMax,rhoMin,rhoMax);
Di_MC->Draw("QGFit4>>QG4_Q",sel,"goff",entries);
sprintf(sel,"eventWeight*(pdgIdJet0!=0 && abs(pdgIdJet0)==5 &&%d<ptJet0 &&ptJet0<%d &&%d<rhoPF &&rhoPF<%d && abs(etaJet0<2))",ptMin,ptMax,rhoMin,rhoMax);
Di_MC->Draw("QGFit4>>QG4_B",sel,"goff",entries);
sprintf(sel,"eventWeight*( ( pdgIdJet0==0 ||(abs(pdgIdJet0)>5 &&pdgIdJet0!=21) ) &&%d<ptJet0 &&ptJet0<%d &&%d<rhoPF &&rhoPF<%d && abs(etaJet0<2))",ptMin,ptMax,rhoMin,rhoMax);
Di_MC->Draw("QGFit4>>QG4_U",sel,"goff",entries);

//Scale
TH1F*all=QG4_Q->Clone("all");
all->Add(QG4_G);
all->Add(QG4_B);
all->Add(QG4_U);
printf("%.2f %.2f %.2f %.2f =%.2f\n",QG4_Q->Integral(),QG4_G->Integral(),QG4_B->Integral(),QG4_U->Integral(),all->Integral());

QG4_Q->Scale(1./all->Integral("width"));
QG4_G->Scale(1./all->Integral("width"));
QG4_B->Scale(1./all->Integral("width"));
QG4_U->Scale(1./all->Integral("width"));

//add to a Stack
THStack *S=new THStack("stack","stack");
QG4_Q->GetXaxis()->SetTitle("Q-G Discr.");
QG4_Q->GetYaxis()->SetTitle("Normalized to Unity");

S->Add(QG4_Q);
S->Add(QG4_G);
S->Add(QG4_B);
S->Add(QG4_U);

TCanvas *c2=new TCanvas();
S->Draw();
//all->SetLineColor(kRed);all->SetLineWidth(2);all->Scale(1./all->Integral("width"));all->Draw("HIST SAME");
QG4_Di->Draw("P SAME");

TLegend *L=new TLegend(0.4,.5,0.6,.7);
L->AddEntry("QG4_Di","QG Data","PF");
L->AddEntry("QG4_Q","QG Quark","F");
L->AddEntry("QG4_G","QG Gluon","F");
L->AddEntry("QG4_B","QG B","F");
L->AddEntry("QG4_U","QG Undefined","F");
L->Draw();
CMSLatex *a=new CMSLatex();
a->DrawSimulation();
TLatex *l=new TLatex();
l->SetNDC();
l->SetTextSize(0.04);
l->SetTextAlign(23);
if(what=='P')l->DrawLatex(0.75,.85,"#gamma Jet");
if(what=='D')l->DrawLatex(0.75,.85,"Di Jet");

c2->RedrawAxis();
if(what=='P')c2->SaveAs("QGFit4_GammaJet.pdf");
if(what=='D')c2->SaveAs("QGFit4_DiJet.pdf");
return 0;
}
