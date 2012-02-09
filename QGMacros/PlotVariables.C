#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
//printf - sprintf - fprintf
#include <stdio.h>
//
#include <stdlib.h>
//=================================================================
void PlotVariables(const char *mcFile, const char *dataFile="",
const char range[]="(50,0,1)",
const char Name[]="(ptDJetReco)",
const char Title[]="Likelihood",
const char outFileName[]="",
const double PtMin=85,
const double PtMax=115,
const double RhoMin=0,
const double RhoMax=20,
const bool DrawGlobal=true,
const char *addSel=""
)
//=================================================================
{
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetOptStat(kFALSE);

bool DrawData=(dataFile[0]!='\0');

TFile *f0=TFile::Open(mcFile);
TFile *f1=NULL;
	if(DrawData)f1=TFile::Open(dataFile);

 TTree *tmc=(TTree*)f0->Get("jetTree");
 TTree *tdata=NULL;
	if(DrawData)tdata=(TTree*)f1->Get("jetTree");

//some options
// utility variables
char str[1023];
char cut[1023];
int nBins;
float xMin,xMax;
sscanf(range,"(%d,%f,%f)",&nBins,&xMin,&xMax);
//common selection
//const char selection[]=" && abs(etaJetReco)<2.4 && passedPhotonID_medium && trackCountingHighEffBJetTagsJetReco<1.7";
 char selection[1023];
	sprintf(selection," && abs(etaJetReco)<2.4 && trackCountingHighEffBJetTagsJetReco<1.7 %s",addSel);
//===============================================GETTING HISTOS===================================

// DiJet MC Quark
 TH1F* LDiMCQ=new TH1F("LDiMCQ","LDiMCQ",nBins,xMin,xMax);
 LDiMCQ->Sumw2();
 sprintf(str,"%s>>LDiMCQ",Name);
 sprintf(cut,"eventWeight*(%lf<ptJetReco && ptJetReco <%lf && abs(pdgIdPart)<4  && %lf<rhoPF && rhoPF<%lf %s ) ",PtMin,PtMax,RhoMin,RhoMax,selection);
 tmc->Draw(str,cut,"goff E");
 //tmc->Draw(str,cut," E");
 LDiMCQ->SetLineColor(kBlack);
 LDiMCQ->SetFillColor(kGray);
 LDiMCQ->GetXaxis()->SetTitle(Title);
 LDiMCQ->Scale(1./LDiMCQ->Integral());

//DiJet MC Gluon
 TH1F* LDiMCG=new TH1F("LDiMCG","LDiMCG",nBins,xMin,xMax);
 LDiMCG->Sumw2();
 sprintf(str,"%s>>LDiMCG",Name);
 sprintf(cut,"eventWeight*(%lf<ptJetReco && ptJetReco <%lf && abs(pdgIdPart)==21 && %lf<rhoPF && rhoPF<%lf %s)",PtMin,PtMax,RhoMin,RhoMax,selection);
 tmc->Draw(str,cut,"goff E");
 //tmc->Draw(str,cut," E");
 LDiMCG->SetLineColor(kRed);
 LDiMCG->SetFillColor(kRed-9);
 LDiMCG->Scale(1./LDiMCG->Integral());

//DiJet MC
 TH1F* LDiMC=new TH1F("LDiMC","LDiMC",nBins,xMin,xMax);
 LDiMC->Sumw2();
 sprintf(str,"%s>>LDiMC",Name);
 sprintf(cut,"(%lf<ptJetReco && ptJetReco <%lf &&%lf<rhoPF && rhoPF<%lf %s)*eventWeight",PtMin,PtMax,RhoMin,RhoMax,selection);
 tmc->Draw(str,cut,"goff E");
 LDiMC->SetLineColor(kGreen+3);
 LDiMC->SetFillColor(kGreen-9);
 LDiMC->Scale(1./LDiMC->Integral());


//========== DATI=============


//Di
 TH1F* LDidata=new TH1F("LDidata","LDidata",nBins,xMin,xMax);
if(DrawData){
 sprintf(str,"%s>>LDidata",Name);
 sprintf(cut,"%lf<ptJetReco && ptJetReco <%lf && %lf<rhoPF && rhoPF<%lf %s",PtMin,PtMax,RhoMin,RhoMax,selection);
 tdata->Draw(str,cut,"goff E");
 //TH1F*LDidata=gDirectory->Get("LDidata")->Clone();
 LDidata->SetMarkerStyle(22);
 LDidata->SetMarkerColor(kMagenta);
 LDidata->SetMarkerSize(1.0);
 LDidata->Sumw2();
 LDidata->Scale(1./LDidata->Integral());
}



//========================
//Confronto MC dati Ph+jet
//========================

float Max=0;
Max=(Max>LDiMCQ->GetMaximum())?Max:LDiMCQ->GetMaximum();
Max=(Max>LDiMCG->GetMaximum())?Max:LDiMCG->GetMaximum();
Max=(Max>LDiMC->GetMaximum())?Max:LDiMC->GetMaximum();
if(DrawData)Max=(Max>LDidata->GetMaximum())?Max:LDidata->GetMaximum();
Max*=1.1;
LDiMCQ->SetMaximum(Max);

TCanvas *c1=new TCanvas("c2");
LDiMCQ->GetYaxis()->SetTitle("Normalized to Unity");
LDiMCQ->Draw("HIST");
LDiMCG->Draw("HIST SAME");
TH1F*tmp=(TH1F*)LDiMCQ->Clone("tmp");
tmp->SetFillColor(0);
tmp->Draw("HIST SAME");
LDiMC->SetFillStyle(3004);
if(DrawGlobal)LDiMC->Draw("HIST SAME");
if(DrawData)LDidata->Draw("P E 0 SAME");

sprintf(str,"%.0lf < P_{T}[GeV/c]< %.0lf     %.0lf<#rho<%.0lf",PtMin,PtMax,RhoMin,RhoMax);
TLegend *L=new TLegend(0.6,0.6,0.8,.89,str);
L->SetFillColor(0);
L->SetBorderSize(0);
if(DrawData)L->AddEntry("LDidata","#gamma+jet data","LP");
//L->AddEntry("LDiMC","#gamma+jet MonteCarlo","F");
if(DrawGlobal)L->AddEntry("LDiMC","MonteCarlo (nat. mixing)","F");
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

if(outFileName[0]!='\0'){
	c1->SaveAs(outFileName);
	//char str[1023];
	sprintf(str,"%s.root",outFileName);
	c1->SaveAs(str);
	}
}
