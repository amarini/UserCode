#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TH1F.h"
//printf - sprintf - fprintf
#include <stdio.h>
//
#include <stdlib.h>
//.h -> only headers
#include "QGLikelihood/QGLikelihoodCalculator.C"

//root Di_mc_88_110.root Di_data_B_88_110.root Ph_mc_88_110.root Ph_data_B_88_110.root
void PlotChiSquare(const char *mcFile, const char *dataFile,
const char range[]="(50,0,1)",
const char Name[]="(ptDJetReco)",
const char Title[]="Likelihood",
const double PtMin=85,
const double PtMax=115
)
{
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetOptStat(kFALSE);

TFile *f0=TFile::Open(mcFile);
TFile *f1=TFile::Open(dataFile);

 TTree *tmc=f0->Get("jetTree");
 TTree *tdata=f1->Get("jetTree");

//some options
// utility variables
char str[1023];
char cut[1023];
int nBins;
float xMin,xMax;
sscanf(range,"(%d,%f,%f)",&nBins,&xMin,&xMax);
//common selection
const char selection[]=" && abs(etaJetReco)<2.4 && passedPhotonID_medium && trackCountingHighEffBJetTagsJetReco<1.7";
//===============================================GETTING HISTOS===================================

// DiJet MC Quark
 TH1F* LDiMCQ=new TH1F("LDiMCQ","LDiMCQ",nBins,xMin,xMax);
 LDiMCQ->Sumw2();
 sprintf(str,"%s>>LDiMCQ",Name);
 sprintf(cut,"eventWeight*(%lf<ptJetReco && ptJetReco <%lf && abs(pdgIdPart)<4) %s ",PtMin,PtMax,selection);
 tmc->Draw(str,cut,"goff E");
 LDiMCQ->SetLineColor(kBlack);
 LDiMCQ->SetFillColor(kGray);
 LDiMCQ->GetXaxis()->SetTitle(Title);
 LDiMCQ->Scale(1./LDiMCQ->Integral());

//DiJet MC Gluon
 TH1F* LDiMCG=new TH1F("LDiMCG","LDiMCG",nBins,xMin,xMax);
 LDiMCG->Sumw2();
 sprintf(str,"%s>>LDiMCG",Name);
 sprintf(cut,"eventWeight*(%lf<ptJetReco && ptJetReco <%lf && abs(pdgIdPart)==21 %s)",PtMin,PtMax,selection);
 tmc->Draw(str,cut,"goff E");
 LDiMCG->SetLineColor(kRed);
 LDiMCG->SetFillColor(kRed-9);
 LDiMCG->Scale(1./LDiMCG->Integral());

//DiJet MC
 TH1F* LDiMC=new TH1F("LDiMC","LDiMC",nBins,xMin,xMax);
 LDiMC->Sumw2();
 sprintf(str,"%s>>LDiMC",Name);
 sprintf(cut,"(%lf<ptJetReco && ptJetReco <%lf %s)*eventWeight",PtMin,PtMax,selection);
 tmc->Draw(str,cut,"goff E");
 LDiMC->SetLineColor(kGreen);
 LDiMC->SetFillColor(kGreen-9);
 LDiMC->Scale(1./LDiMC->Integral());


//========== DATI=============


//Di
 TH1F* LDidata=new TH1F("LDidata","LDidata",nBins,xMin,xMax);
 sprintf(str,"%s>>LDidata",Name);
 sprintf(cut,"%lf<ptJetReco && ptJetReco <%lf %s",PtMin,PtMax,selection);
 tdata->Draw(str,cut,"goff E");
 //TH1F*LDidata=gDirectory->Get("LDidata")->Clone();
 LDidata->SetMarkerStyle(22);
 LDidata->SetMarkerColor(kMagenta);
 LDidata->SetMarkerSize(1.0);
 LDidata->Sumw2();
 LDidata->Scale(1./LDidata->Integral());



//=========================================================================DRAW=============================

TGraph *G=new TGraph();
int count=0;

TCanvas *c1= new TCanvas;
//I want to Draw the minimum chi square and the ChiSquare Dependence
c1->Divide(2);
c1->cd(1);
double MinChi2=1e10;//infty

for(double frac=1.0;frac>=0;frac-=0.005){
sprintf(str,"newq%.2lf",frac);
TH1F*hq=(TH1F*)LDiMCQ->Clone(str);

sprintf(str,"newg%.2lf",frac);
TH1F*hg=(TH1F*)LDiMCG->Clone(str);

hq->Scale(frac);
hg->Scale(1.0-frac);
hq->Add(hg);
double Chi2=LDidata->Chi2Test(hq,"WW CHI2/NDF");
cout<<"frac quark "<<frac<<" - Chi2 "<< Chi2<<endl;
G->SetPoint(count,frac,Chi2);count++;
hq->SetLineColor(kRed);
hq->SetFillColor(kOrange);
if(Chi2<MinChi2)hq->Draw("HIST 0");
if(Chi2<MinChi2)LDidata->Draw("HIST P E 0 SAME");
sprintf(str,"frac quark = %.2lf",frac);
TLegend *L=new TLegend(0.6,0.6,0.8,.89,str);
sprintf(str,"newq%.2lf",frac);
L->AddEntry(str,"MC","F");
L->AddEntry("LDidata","Dijet data","PLF");
if(Chi2<MinChi2)L->Draw();
if(Chi2<MinChi2)MinChi2=Chi2;
}
G->GetYaxis()->SetRangeUser(MinChi2/2.0,MinChi2*10);
//
//
c1->cd(2);
G->GetXaxis()->SetTitle("Quark Fraction");
G->GetYaxis()->SetTitleOffset(1.2);
G->GetYaxis()->SetTitle("#chi^{2}/NDF");
G->Draw("ALP");

/**************PhotonData********/
//for(int k=1;k<=LPhMC->GetNbinsX();k++){LPhMC->SetBinError(k,0);}
//for(int k=1;k<=LPhMCQ->GetNbinsX();k++){LPhMCQ->SetBinError(k,0);}
//for(int k=1;k<=LPhMCG->GetNbinsX();k++){LPhMCG->SetBinError(k,0);}
//
//TGraph *G1=new TGraph();
//
//int count=0;
//
//c1=new TCanvas();
//c1->Divide(2,2);
//i=1;
//
//for(double frac=0.90;frac>=0.65;frac-=0.02){
//c1->cd(i);
//if((0.815<frac&&frac<0.825)||(.855<frac&&frac<0.865)||(.775<frac&&frac<.785)){i++;}
//sprintf(str,"newq%.2lf",frac);
//TH1F*hq=(TH1F*)LPhMCQ->Clone(str);
//
//sprintf(str,"newg%.2lf",frac);
//TH1F*hg=(TH1F*)LPhMCG->Clone(str);
//
//hq->Scale(frac);
//hg->Scale(1.0-frac);
//hq->Add(hg);
//double Chi2=LPhdata->Chi2Test(hq,"WW CHI2/NDF");
//cout<<"frac quark "<<frac<<" - Chi2 "<< Chi2<<endl;
//G1->SetPoint(count,frac,Chi2);count++;
//hq->SetLineColor(kRed);
//hq->SetFillColor(kOrange);
//hq->Draw("HIST 0");
//LPhdata->Draw("HIST P E 0 SAME");
//sprintf(str,"frac quark = %.2lf",frac);
//TLegend *L=new TLegend(0.6,0.6,0.8,.89,str);
//sprintf(str,"newq%.2lf",frac);
//L->AddEntry(str,"MC","F");
//L->AddEntry("LPhdata","#gamma+jet data","PLF");
//L->Draw();
//}
//
//c1->cd(4);
//G1->GetXaxis()->SetTitle("Quark Fraction");
//G1->GetYaxis()->SetTitleOffset(1.2);
//G1->GetYaxis()->SetTitle("#chi^{2}/NDF");
//G1->Draw("ALP");

//G1->SetPoint(count,(frac-0.42)/(0.87-0.42),Chi2);count++;

//========================
//Confronto MC dati Ph+jet
//========================
TCanvas *c1=new TCanvas("c2");
LDiMCQ->GetYaxis()->SetTitle("Normalized to Unity");
LDiMCQ->Draw("HIST");
LDiMCG->Draw("HIST SAME");
TH1F*tmp=LDiMCQ->Clone("tmp");
tmp->SetFillColor(0);
tmp->Draw("HIST SAME");
LDiMC->SetFillStyle(3004);
LDiMC->Draw("HIST SAME");
LDidata->Draw("P E 0 SAME");

sprintf(str,"%.0lf < P_{T}[GeV/c]< %.0lf",PtMin,PtMax);
TLegend *L=new TLegend(0.6,0.6,0.8,.89,str);
L->SetFillColor(0);
L->SetBorderSize(0);
L->AddEntry("LDidata","#gamma+jet data","LP");
L->AddEntry("LDiMC","#gamma+jet MonteCarlo","F");
L->AddEntry("LDiMCQ","Quark MonteCarlo","F");
L->AddEntry("LDiMCG","Gluon MonteCarlo","F");
L->Draw();
c1->RedrawAxis();


TLatex *lat=new TLatex();
lat->SetNDC();//normalized coordinate syst ref to pad
lat->SetTextSize(0.04);
lat->SetTextAlign(23);
lat->DrawLatex(0.25,.88,"CMS Preliminary");
lat->DrawLatex(0.25,.82,"#sqrt{s} = 7 TeV");


}
