//Author: Andrea Carlo Marini
//email: andrea.carlo.marini@cern.ch
//date: 24/01/2011
#include <stdio.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TROOT.h"

void mkCorrelationPlot(const char * filename)
{
gROOT->SetBatch();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(kFALSE);
#include "PtBins.h"
TFile *f=TFile::Open(filename);
TTree *t=(TTree*)f->Get("demo/t");
char targetfile[255];
char cut[255];
TH2F *h;
TFile *f1;
for(int i=0;i<NBins-1;i++)
	{
	sprintf(targetfile,"CorrelationPlot/Correlation_%03.0lf_%03.0lf.root",PtBins[i],PtBins[i+1]);

	sprintf(cut,"%10.0lf<jtpt && jtpt<=%10.0lf && n>1 && -2.0<jteta && jteta<2.0",PtBins[i],PtBins[i+1]);

	f1=new TFile(targetfile,"RECREATE");
	f1->cd();

	t->Draw("PtD:ncharged>>h1",cut,"box");
	h=(TH2F*)gDirectory->Get("h1");
	h->SetName("PtD vs ncharged");
	h->SetDirectory(f1->GetDirectory(""));
	h->Write();

	t->Draw("PtD:nneutral>>h1",cut,"box");
	h=(TH2F*)gDirectory->Get("h1");
	h->SetName("PtD vs nneutral");
	h->SetDirectory(f1->GetDirectory(""));
	h->Write();

	t->Draw("PtD:r>>h1",cut,"box");
	h=(TH2F*)gDirectory->Get("h1");
	h->SetName("PtD vs rRMS");
	h->SetDirectory(f1->GetDirectory(""));
	h->Write();
	
	t->Draw("ncharged:nneutral>>h1",cut,"box");
	h=(TH2F*)gDirectory->Get("h1");
	h->SetName("Ncharged vs Nneutral");
	h->SetDirectory(f1->GetDirectory(""));
	h->Write();

	t->Draw("ncharged:r>>h1",cut,"box");
	h=(TH2F*)gDirectory->Get("h1");
	h->SetName("ncharged vs rRMS");
	h->SetDirectory(f1->GetDirectory(""));
	h->Write();

	t->Draw("nneutral:r>>h1",cut,"box");
	h=(TH2F*)gDirectory->Get("h1");
	h->SetName("Nneutral vs rRMS");
	h->SetDirectory(f1->GetDirectory(""));
	h->Write();

	sprintf(cut,"pdgid==21 && %lf<jtpt && jtpt<=%lf && n>1 && -2.0<jteta && jteta<2.0",PtBins[i],PtBins[i+1]);


	sprintf(cut,"pdgid!=21 && %lf<jtpt && jtpt<=%lf && n>1 && -2.0<jteta && jteta<2.0",PtBins[i],PtBins[i+1]);

	
	f1->Close();
	}
f->Close();
}
