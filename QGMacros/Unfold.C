#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include <stdlib.h>
#include <stdio.h>
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "RooUnfold/RooUnfoldResponse.h"
#ifdef __CINT__
gSystem->Load("libRooUnfold.so");
#endif

TH1F* MergeHistos(TH1F* a,TH1F*b,const char*name="merged")
{
if(a->GetNbinsX()!=b->GetNbinsX())return NULL;
float aMin=a->GetBinLowEdge(1);
float aMax=a->GetBinLowEdge(a->GetNbinsX()+1);
float bMin=b->GetBinLowEdge(1);
float bMax=b->GetBinLowEdge(b->GetNbinsX()+1);
//fprintf(stderr,"%d: %.5f %.5f\n",a->GetNbinsX(),aMin,aMax);
TH1F*R=new TH1F(name,name,(a->GetNbinsX())*2,aMin,aMax+(aMax-aMin));
R->Sumw2();
for(int i=1;i<=a->GetNbinsX();i++)
	{
	R->SetBinContent(i,a->GetBinContent(i));	
	R->SetBinError(i,a->GetBinError(i));	
	}
for(int i=1;i<=b->GetNbinsX();i++)
	{
	R->SetBinContent(i+a->GetNbinsX(),b->GetBinContent(i));	
	R->SetBinError(i+a->GetNbinsX(),b->GetBinError(i));	
	}
return R;
}

int DivideHisto(TH1F* a,TH1F& *b TH1F& *c)
{
if(  (a->GenBinsX()&1)!=0 ) return 1; //not event number of bins
char name[1023];
sprintf(name,"%s_1",a->GetName());
b=new TH1F(name,name,a->GetNBinsX()/2,a->GetBinLowEdge(1),a->GetBinLowEdge(a->GetNbinsX()+1));
sprintf(name,"%s_2",a->GetName());
c=new TH1F(name,name,a->GetNBinsX()/2,a->GetBinLowEdge(1),a->GetBinLowEdge(a->GetNbinsX()+1));
//TODO	
return 0;
}


int Unfold(const char*fileName1,
	  const char*fileName2,
  	  const char *varName="ptDJet0",
	  const char *range="(100,0,1)",
	  float PtMin=100,float PtMax=110,
	  float RhoMin=12,float RhoMax=13,
	  float alpha1=0.85,float alpha2=0.55,
	  const char*treeName="tree_passedEvents")
{

 //Some stuff
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
//
TFile *f1=TFile::Open(fileName1);
TFile *f2=TFile::Open(fileName2);
if( (f1==NULL) || (f2==NULL)){fprintf(stderr,"FILES DOES NOT EXIST!\n");return 1;}

//Create the Unfold Matrix - merge of the histograms
//use alpha1 alpha2 
	int nBins=100;
	float xMin,xMax;
	sscanf(range,"(%d,%f,%f)",&nBins,&xMin,&xMax);
	TH1F *aux=new TH1F("aux","aux",nBins*2,xMin,xMax+(xMax-xMin));
	RooUnfoldResponse *R=new RooUnfoldResponse(aux,aux,"Unfold");
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(i)      ,aux->GetBinCenter(i)     ,alpha1); //RECO - GEN (PhDi - QG)
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(nBins+i),aux->GetBinCenter(nBins+i),alpha2);
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(nBins+i),aux->GetBinCenter(i)      ,1-alpha1);//TODO: tocheck
	for(int i=1;i<=nBins;i++)R->Fill(aux->GetBinCenter(i)      ,aux->GetBinCenter(nBins+i),1-alpha2);
	TCanvas *c1=new TCanvas();
	((TH2F*)R->Hresponse())->Draw("BOX");
//Getting the histograms
	TTree *t1=(TTree*)f1->Get(treeName);
	TTree *t2=(TTree*)f2->Get(treeName);
	char name[1023];
	char selection[1023];
	TH1F *h1=new TH1F("var1","var1",nBins,xMin,xMax);h1->Sumw2();
	TH1F *h2=new TH1F("var2","var2",nBins,xMin,xMax);h2->Sumw2();
	sprintf(name,"%s>>var1",varName);
	//sprintf(selection," %f < ptJet0 && ptJet0<%f && %f<rho && rho <%f ",PtMin,PtMax,RhoMin,RhoMax);
	sprintf(selection," %f < ptJet0 && ptJet0<%f ",PtMin,PtMax,RhoMin,RhoMax);
	t1->Draw(name,selection,"goff");h1->Scale(1./h1->Integral());

	sprintf(name,"%s>>var2",varName);
	sprintf(selection," %f < ptJet0 && ptJet0<%f && passedPhotonID && !btagged ",PtMin,PtMax,RhoMin,RhoMax);
	t2->Draw(name,selection,"goff");h2->Scale(1./h2->Integral());
//Unfold the distributions
	TH1F* H=MergeHistos(h1,h2,"merged");
	RooUnfoldSvd *Unfold=new RooUnfoldSvd(R,H,nBins);
//	RooUnfoldInvert *Unfold=new RooUnfoldInvert(R,H);
//Print the output
	TH1F *hq=new TH1F("quark","quark",nBins,xMin,xMax);hq->Sumw2();
	sprintf(name,"%s>>quark",varName);
	sprintf(selection," %f < ptJet0 && ptJet0<%f && abs(pdgIdPartJet0) <4 ",PtMin,PtMax,RhoMin,RhoMax);
	t1->Draw(name,selection,"goff");hq->Scale(1./hq->Integral());
	TH1F *hg=new TH1F("gluon","gluon",nBins,xMin,xMax);hg->Sumw2();
	sprintf(name,"%s>>gluon",varName);
	sprintf(selection," %f < ptJet0 && ptJet0<%f && abs(pdgIdPartJet0) ==21 ",PtMin,PtMax,RhoMin,RhoMax);
	t1->Draw(name,selection,"goff");hg->Scale(1./hg->Integral());

	TCanvas *c=new TCanvas();
	sprintf(selection," %f < ptJet0 && ptJet0<%f ",PtMin,PtMax,RhoMin,RhoMax);
	TH1F* unf=Unfold->Hreco(3);
	unf->SetLineColor(kBlue+2);

	TH1F*allTrue=MergeHistos(hg,hq,"true-g-q");
	allTrue->SetLineColor(kRed+2);

	allTrue->Draw("");
	unf->Draw("HIST SAME"); //unfold
	H->SetMarkerColor(kMagenta+2);
	H->SetMarkerStyle(20);
	H->Draw("P SAME");
		
	TLegend *L=new TLegend(0.45,0.75,0.55,.89);
	L->AddEntry("merged","merged");
	L->AddEntry("true-g-q","g-q");
	L->AddEntry("Unfold","Unfolded");
	L->Draw();
return 0;
}

