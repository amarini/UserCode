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


int ComputeMixture(const char*fileName1,
	  const char*fileName2,
	  float PtMin=100,float PtMax=110,
	  float RhoMin=12,float RhoMax=13,
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

//Getting the histograms
	TTree *t1=(TTree*)f1->Get(treeName);
	TTree *t2=(TTree*)f2->Get(treeName);
	const char Jet[]="Jet0";
	char selection[1023];
	float q,g,o;
	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(pdgIdPartJet0)<4)",PtMin,PtMax,RhoMin,RhoMax); //dijet - file 1
	TH1F *hq=new TH1F("quark1","quark1",100,0,100);
		t1->Draw("abs(pdgIdPartJet0)>>quark1",selection,"goff");
		q=hq->Integral();
	//q=t1->GetEntries(selection);
	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(pdgIdPartJet0)==21)",PtMin,PtMax,RhoMin,RhoMax); //gluon
	TH1F *hg=new TH1F("gluon1","gluon1",100,0,100);
		t1->Draw("abs(pdgIdPartJet0)>>gluon1",selection,"goff");
		g=hg->Integral();
	//g=t1->GetEntries(selection);	
	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(pdgIdPartJet0)!=21 && abs(pdgIdPartJet0)>=4)",PtMin,PtMax,RhoMin,RhoMax); //gluon
	TH1F *ho=new TH1F("other1","other1",100,0,100);
		t1->Draw("abs(pdgIdPartJet0)>>other1",selection,"goff");
		o=ho->Integral();
	//o=t1->GetEntries(selection);
	
	fprintf(stderr,"DIJET q/q+g=%.3f q/q+g+o=%.3f o/q+g+o=%.3f\n",q/(q+g),q/(q+g+o),o/(q+g+o))	;

	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF <%f&& passedID_FULL && !btagged&& abs(pdgIdPartJet0)<4)",PtMin,PtMax,RhoMin,RhoMax);
	 	hq=new TH1F("quark2","quark2",100,0,100);
		t2->Draw("abs(pdgIdPartJet0)>>quark2",selection,"goff");
		q=hq->Integral();

	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF <%f&& passedID_FULL && !btagged && abs(pdgIdPartJet0)==21)",PtMin,PtMax,RhoMin,RhoMax);
		hg=new TH1F("gluon2","gluon2",100,0,100);
		t2->Draw("abs(pdgIdPartJet0)>>gluon2",selection,"goff");
		g=hg->Integral();

	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF <%f&& passedID_FULL && !btagged && abs(pdgIdPartJet0)!=21 && abs(pdgIdPartJet0)>=4)",PtMin,PtMax,RhoMin,RhoMax);
	TH1F *ho=new TH1F("other2","other2",100,0,100);
		t2->Draw("abs(pdgIdPartJet0)>>other2",selection,"goff");
		o=ho->Integral();

	fprintf(stderr,"PHJET q/q+g=%.3f q/q+g+o=%.3f o/q+g+o=%.3f\n",q/(q+g),q/(q+g+o),o/(q+g+o))	;
TCanvas *c1=new TCanvas();
c1->Divide(2);
c1->SetLogx();
	c1->cd(1);
	t1->Draw("rhoPF:ptJet0>>dijet1","","BOX");
	t2->Draw("rhoPF:ptJet0>>phjet1","passedID_FULL && !btagged","BOX SAME");
	TH1F*ph=gDirectory->Get("phjet1");
	ph->SetFillColor(kRed);
	ph->SetLineColor(0);
	c1->cd(2);
	t1->Draw("rhoPF:ptJet0>>dijet2","eventWeight","BOX");
	t2->Draw("rhoPF:ptJet0>>phjet2","eventWeight*(passedID_FULL && !btagged)","BOX SAME");
	TH1F*ph2=gDirectory->Get("phjet2");
	ph2->SetFillColor(kRed);
	ph2->SetLineColor(0);
return 0;
}

