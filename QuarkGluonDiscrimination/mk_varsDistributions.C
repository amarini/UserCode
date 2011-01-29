//Author: Andrea Carlo Marini
//email: andrea.carlo.marini@cern.ch
//date: 24/01/2011
#include <stdio.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TROOT.h"

void mkvarsdistributions(const char * filename)
{
gROOT->SetBatch();
#include "PtBins.h"
TFile *f=TFile::Open(filename);
TTree *t=(TTree*)f->Get("demo/t");
char targetfile[255];
char cut[255];
TH1F *h;
TFile *f1;
for(int i=0;i<NBins-1;i++)
	{
	sprintf(targetfile,"VarDistribution_%03.0lf_%03.0lf.root",PtBins[i],PtBins[i+1]);

	sprintf(cut,"%10.0lf<jtpt && jtpt<=%10.0lf && n>1 && -2.0<jteta && jteta<2.0",PtBins[i],PtBins[i+1]);

	f1=new TFile(targetfile,"RECREATE");
	f1->cd();

	t->Draw("ncharged>>h1(100,0,100)",cut);
	h=(TH1F*)gDirectory->Get("h1");
	h->SetName("ncharged");
	h->SetDirectory(f1->GetDirectory(""));
	h->Write();

	t->Draw("nneutral>>h2(100,0,100)",cut);
	h=(TH1F*)gDirectory->Get("h2");
	h->SetName("nneutral");
	h->Write();
	
	t->Draw("PtD>>h3(100,0,1)",cut);
	h=(TH1F*)gDirectory->Get("h3");
	h->SetName("PtD");
	h->Write();
	
	t->Draw("r>>h4(1000,0,0.1)",cut);
	h=(TH1F*)gDirectory->Get("h4");
	h->SetName("rRMS");
	h->Write();

	sprintf(cut,"pdgid==21 && %lf<jtpt && jtpt<=%lf && n>1 && -2.0<jteta && jteta<2.0",PtBins[i],PtBins[i+1]);

	t->Draw("ncharged>>h5(100,0,100)",cut);
	h=(TH1F*)gDirectory->Get("h5");
	h->SetName("ncharged_gluon");
	h->Write();

	t->Draw("nneutral>>h6(100,0,100)",cut);
	h=(TH1F*)gDirectory->Get("h6");
	h->SetName("nneutral_gluon");
	h->Write();
	
	t->Draw("PtD>>h7(100,0,1)",cut);
	h=(TH1F*)gDirectory->Get("h7");
	h->SetName("PtD_gluon");
	h->Write();
	
	t->Draw("r>>h8(1000,0,0.1)",cut);
	h=(TH1F*)gDirectory->Get("h8");
	h->SetName("rRMS_gluon");
	h->Write();

	sprintf(cut,"pdgid!=21 && %lf<jtpt && jtpt<=%lf && n>1 && -2.0<jteta && jteta<2.0",PtBins[i],PtBins[i+1]);

	t->Draw("ncharged>>h11(100,0,100)",cut);
	h=(TH1F*)gDirectory->Get("h11");
	h->SetName("ncharged_quark");
	h->Write();

	t->Draw("nneutral>>h12(100,0,100)",cut);
	h=(TH1F*)gDirectory->Get("h12");
	h->SetName("nneutral_quark");
	h->Write();
	
	t->Draw("PtD>>h13(100,0,1)",cut);
	h=(TH1F*)gDirectory->Get("h13");
	h->SetName("PtD_quark");
	h->Write();
	
	t->Draw("r>>h14(1000,0,0.1)",cut);
	h=(TH1F*)gDirectory->Get("h14");
	h->SetName("rRMS_quark");
	h->Write();
	
	f1->Close();
	}
f->Close();
}
