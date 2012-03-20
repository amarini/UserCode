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
#include <iostream>
using namespace std;
int ComputeMixture(const char*fileName1,
	  float PtMin=100,float PtMax=110,
	  float RhoMin=12,float RhoMax=13,
	  const char*treeName="tree_passedEvents",
	  double *x=NULL,double *ex=NULL,
	  double *y=NULL,double *ey=NULL)
{

 //Some stuff
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
//
	TChain *t1=new TChain(treeName);
	cout << t1->Add(fileName1)<<" added "<<endl;

	char selection[1023];
	float q,g,o;
	float eq,eg,eo;
	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(pdgIdJet0)<=4)",PtMin,PtMax,RhoMin,RhoMax); //dijet - file 1

	TH1F *hq=new TH1F("quark1","quark1",100,0,100);
		hq->Sumw2();
		t1->Draw("abs(pdgIdJet0)>>quark1",selection,"goff");// I'm interested only in the integral
		q=hq->IntegralAndError(1,100,eq);
		fprintf(stderr,"selection quark ==%s==\nq=%f\n",selection,q);

	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(pdgIdJet0)==21)",PtMin,PtMax,RhoMin,RhoMax); //gluon
	TH1F *hg=new TH1F("gluon1","gluon1",100,0,100);
		hg->Sumw2();
		t1->Draw("abs(pdgIdJet0)>>gluon1",selection,"goff");
		g=hg->IntegralAndError(1,100,eg);
		fprintf(stderr,"selection gluon ==%s==\ng=%f\n",selection,g);

	sprintf(selection,"eventWeight*( %f < ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(pdgIdJet0)!=21 && abs(pdgIdJet0)>4 )",PtMin,PtMax,RhoMin,RhoMax); //gluon
	TH1F *ho=new TH1F("other1","other1",100,0,100);
		ho->Sumw2();
		t1->Draw("abs(pdgIdJet0)>>other1",selection,"goff");
		o=ho->IntegralAndError(1,100,eo);
		fprintf(stderr,"selection others ==%s==\no=%f\n",selection,o);
	
	fprintf(stderr,"DIJET q/q+g=%.3f q/q+g+o=%.3f o/q+g+o=%.3f\n",q/(q+g),q/(q+g+o),o/(q+g+o));
	if(x!=NULL) *x=q/(q+g+o);
	if(x!=NULL && ex!=NULL) *ex=TMath::Sqrt( (eq*eq+eg*eg+eo*eo)/((q+g+o)*(q+g+o)) + eq*eq/(q*q) )*(q/(q+g+o));
	if(y!=NULL) *y=g/(q+g+o);
	if(y!=NULL && ey!=NULL) *ey=TMath::Sqrt( (eq*eq+eg*eg+eo*eo)/((q+g+o)*(q+g+o)) + eg*eg/(g*g) )*(g/(q+g+o));

return 0;
}

