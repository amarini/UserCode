//Author: Andrea Carlo Marini
//email: andrea.carlo.marini@cern.ch
//date: 24/01/2011
//this file uses the VarsDistributions computed from bla bla in order to 
//make the likelihood distribution for a test sample
#include <stdio.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "QuarkGluonDiscrimination.C"//link a tutti gli header e codice
void mklikelihood(const char*filename)
{
gROOT->SetBatch();
//const int NBins=13;
//Double_t PtBins[]={15,30,50,80,100,120,150,170,200,210,300,470,600};
#include "PtBins.h"
TFile *f=TFile::Open(filename);
TTree *t=(TTree*)f->Get("demo/t");
char targetfile[255];
TH1F *h,*hq,*hg;
TFile *f1;
Discrimination A;
Double_t PtD,rRMS,jtpt,jteta;
Int_t ncharged,nneutral,pdgid;
t->SetBranchAddress("PtD",&PtD);
t->SetBranchAddress("r",&rRMS);
t->SetBranchAddress("ncharged",&ncharged);
t->SetBranchAddress("nneutral",&nneutral);
t->SetBranchAddress("pdgid",&pdgid);
t->SetBranchAddress("jtpt",&jtpt);
t->SetBranchAddress("jteta",&jteta);

Double_t R;
for(int i=0;i<NBins-1;i++)
	{
	sprintf(targetfile,"LikelihoodDistribution_%03.0lf_%03.0lf.root",PtBins[i],PtBins[i+1]);
	f1=new TFile(targetfile,"RECREATE");
	f1->cd();
	h=new TH1F("Likelihood","Likelihood",100,0,1.);
	h->SetLineColor(kGreen);
	hq=new TH1F("Likelihood_quark","Likelihood_quark",100,0,1.);
	hg=new TH1F("Likelihood_gluon","Likelihood_gluon",100,0,1.);
	hg->SetLineColor(kRed);
	
	h->SetDirectory(f1->GetDirectory(""));
	hq->SetDirectory(f1->GetDirectory(""));
	hg->SetDirectory(f1->GetDirectory(""));
	
	//Loop on the tree
	for(Int_t j =0 ; j<t->GetEntries();j++)	
		{
		t->GetEntry(j);
		//selection of event
		if((PtBins[i+1]<jtpt)||(jtpt<PtBins[i])) continue;//pt binnig selection
		if((jteta<-2.0) || (jteta>2.0))continue; //selection on eta 
		if(nneutral+ncharged<=1)continue; //selection on n
		R=A.IndipendentLikelihood(jtpt,ncharged,nneutral,PtD,rRMS);
		h->Fill(R);
		if(pdgid==21)hg->Fill(R);
		else hq->Fill(R);
		}
	h->Write();
	hq->Write();
	hg->Write();
	f1->Write();
	f1->Close();
	}
}
