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
//#include "RooUnfold/RooUnfoldResponse.h"
//#ifdef __CINT__
//gSystem->Load("libRooUnfold.so");
//#endif
#include <iostream>
using namespace std;
int ComputeMixture(const char*fileName1,
	  const char*fileName2,
	  float PtMin=100,float PtMax=110,
	  float RhoMin=12,float RhoMax=13,
	  const char jet1[]="Jet0",
	  const char jet2[]="Jet0",
	  const char*treeName="tree_passedEvents")
{

 //Some stuff
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
//
//TFile *f1=TFile::Open(fileName1);
//TFile *f2=TFile::Open(fileName2);
TChain *t1=new TChain(treeName);
TChain *t2=new TChain(treeName);
cout << t1->Add(fileName1)<<" added to DIJET"<<endl;
cout << t2->Add(fileName2)<<" added to PHJET"<<endl;

//Getting the histograms
	//TTree *t1=(TTree*)f1->Get(treeName);
	//TTree *t2=(TTree*)f2->Get(treeName);
	char selection[1023];
	float q,g,o;
	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF<%f && abs(pdgId%s)<4)",PtMin,jet1,jet1,PtMax,RhoMin,RhoMax,jet1); //dijet - file 1
	TH1F *hq=new TH1F("quark1","quark1",100,0,100);
		t1->Draw("abs(pdgIdJet0)>>quark1",selection,"goff");// I'm interested only in the integral
		q=hq->Integral();
		fprintf(stderr,"selection quark ==%s==\nq=%f\n",selection,q);
	//q=t1->GetEntries(selection);
	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF<%f && abs(pdgId%s)==21)",PtMin,jet1,jet1,PtMax,RhoMin,RhoMax,jet1); //gluon
	TH1F *hg=new TH1F("gluon1","gluon1",100,0,100);
		t1->Draw("abs(pdgIdJet0)>>gluon1",selection,"goff");
		g=hg->Integral();
		fprintf(stderr,"selection gluon ==%s==\ng=%f\n",selection,g);
	//g=t1->GetEntries(selection);	
	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF<%f && abs(pdgId%s)!=21 && abs(pdgId%s)>=4 )",PtMin,jet1,jet1,PtMax,RhoMin,RhoMax,jet1,jet1); //gluon
	TH1F *ho=new TH1F("other1","other1",100,0,100);
		t1->Draw("abs(pdgIdJet0)>>other1",selection,"goff");
		o=ho->Integral();
		fprintf(stderr,"selection others ==%s==\no=%f\n",selection,o);
	//o=t1->GetEntries(selection);
	
	fprintf(stderr,"DIJET q/q+g=%.3f q/q+g+o=%.3f o/q+g+o=%.3f\n",q/(q+g),q/(q+g+o),o/(q+g+o))	;

	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF <%f&& abs(pdgId%s)<4)",PtMin,jet2,jet2,PtMax,RhoMin,RhoMax,jet2);
	 	hq=new TH1F("quark2","quark2",100,0,100);
		t2->Draw("abs(pdgIdJet0)>>quark2",selection,"goff");
		q=hq->Integral();
		fprintf(stderr,"selection quark ==%s==\nq=%f\n",selection,q);

	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF <%f && abs(pdgId%s)==21)",PtMin,jet2,jet2,PtMax,RhoMin,RhoMax,jet2);
		hg=new TH1F("gluon2","gluon2",100,0,100);
		t2->Draw("abs(pdgIdJet0)>>gluon2",selection,"goff");
		g=hg->Integral();
		fprintf(stderr,"selection gluon ==%s==\ng=%f\n",selection,g);

	sprintf(selection,"eventWeight*( %f < pt%s && pt%s<%f && %f<rhoPF && rhoPF <%f && abs(pdgId%s)!=21 && abs(pdgId%s)>=4)",PtMin,jet2,jet2,PtMax,RhoMin,RhoMax,jet2,jet2);
	TH1F *ho=new TH1F("other2","other2",100,0,100);
		t2->Draw("abs(pdgIdJet0)>>other2",selection,"goff");
		o=ho->Integral();
		fprintf(stderr,"selection others ==%s==\no=%f\n",selection,o);

	fprintf(stderr,"PHJET q/q+g=%.3f q/q+g+o=%.3f o/q+g+o=%.3f\n",q/(q+g),q/(q+g+o),o/(q+g+o))	;
TCanvas *c1=new TCanvas();
c1->Divide(2);
c1->SetLogx();
	c1->cd(1);
	t1->Draw("rhoPF:ptJet0>>dijet1","","BOX");
	t2->Draw("rhoPF:ptJet0>>phjet1","","BOX SAME");
	TH1F*ph=gDirectory->Get("phjet1");
	ph->SetFillColor(kRed);
	ph->SetLineColor(0);
	c1->cd(2);
	t1->Draw("rhoPF:ptJet0>>dijet2","eventWeight","BOX");
	t2->Draw("rhoPF:ptJet0>>phjet2","eventWeight","BOX SAME");
	TH1F*ph2=gDirectory->Get("phjet2");
	ph2->SetFillColor(kRed);
	ph2->SetLineColor(0);
return 0;
}

