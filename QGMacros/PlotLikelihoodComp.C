
#include "PtBins.h"
#include <stdio.h>

int PlotLikelihoodComp()
{
gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetOptStat(kFALSE);

TFile *f=TFile::Open("/Users/andreamarini/Documents/QGDiscriminator/2ndLevel/pandolf/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2_TREE.root");
TTree *t=(TTree*)f->Get("tree_passedEvents");

double PtBins[100];
getBins(PtBins,20,15,1000,true);

char name[1023];
char sel[1023];

for(int PtBin=0;PtBin<20;PtBin++)
	{
	TCanvas *c1=new TCanvas();
	TH1F* lq=new TH1F("LQuark","LQuark",20,0,1.0001);lq->Sumw2();
	TH1F* lg=new TH1F("LGluon","LGluon",20,0,1.0001);lg->Sumw2();
	TH1F* fq=new TH1F("FQuark","FQuark",20,0,1.0001);fq->Sumw2();
	TH1F* fg=new TH1F("FGluon","FGluon",20,0,1.0001);fg->Sumw2();

	sprintf(sel,"eventWeight*(%f<ptJet0 && ptJet0<%f && abs(pdgIdPartJet0)<5 && rhoPF>4 && rhoPF<17 && abs(etaJet0)<2.0)",PtBins[PtBin],PtBins[PtBin+1]);
	t->Draw("Likelihood>>LQuark",sel);		
	t->Draw("LikelihoodFit>>FQuark",sel);		
	sprintf(sel,"eventWeight*(%f<ptJet0 && ptJet0<%f && abs(pdgIdPartJet0)==21 && rhoPF>4 && rhoPF<17 && abs(etaJet0)<2.0)",PtBins[PtBin],PtBins[PtBin+1]);
	t->Draw("Likelihood>>LGluon",sel);		
	t->Draw("LikelihoodFit>>FGluon",sel);
	//Normalize	
	lq->Scale(1./lq->Integral());
	lg->Scale(1./lg->Integral());
	fq->Scale(1./fq->Integral());
	fg->Scale(1./fg->Integral());
	//SetProperties
	double Max=lq->GetMaximum();
	Max=(Max>lg->GetMaximum())?Max:lg->GetMaximum();
	Max=(Max>fg->GetMaximum())?Max:fg->GetMaximum();
	Max=(Max>fq->GetMaximum())?Max:fq->GetMaximum();
	Max*=1.2;
	lg->SetMaximum(Max);
	lg->GetXaxis()->SetTitle("Likelihood");
		
	lg->SetLineColor(kRed+2);lg->SetFillColor(kRed+2);lg->SetFillStyle(3004);lg->SetLineWidth(2);
	lq->SetLineColor(kBlue+2);lq->SetFillColor(kBlue+2);lq->SetFillStyle(3005);lq->SetLineWidth(2);
	fg->SetMarkerColor(kRed);fg->SetMarkerStyle(20);fg->SetMarkerSize(0.8);fg->SetLineColor(kRed);
	fq->SetMarkerColor(kBlue);fq->SetMarkerStyle(20);fq->SetMarkerSize(0.8);fq->SetLineColor(kBlue);
	
	lg->Draw("AXIS");
	lg->Draw("AXIS X+ Y+ SAME");
	lg->Draw("HIST SAME");
	lq->Draw("HIST SAME");
	fq->Draw("P SAME");
	fg->Draw("P SAME");

	sprintf(name,"%.0f < P_{T} <%.0f",PtBins[PtBin],PtBins[PtBin+1]);
	TLegend *L=new TLegend(0.5-0.15,0.7,0.5+0.15,.89,name);	
	L->AddEntry("LQuark","Quark Likelihood","F");
	L->AddEntry("LGluon","Gluon Likelihood","F");
	L->AddEntry("FQuark","Quark Fit","FP");
	L->AddEntry("FGluon","Gluon Fit","FP");
	L->Draw();
		
	sprintf(name,"Likelihood_pt%.0f_%.0f.pdf",PtBins[PtBin],PtBins[PtBin+1]);
	c1->SaveAs(name);	
	delete lq;
	delete lg;
	delete fq;
	delete fg;
	}//PtBins	
}
