#include <stdio.h>
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
//#include "QuarkGluonDiscrimination.h"
void Plot()
{
gROOT->SetBatch();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(kFALSE);
#include "PtBins.h"
char  targetfile[255];
TFile *f;
TH1F *h,*hq,*hg,*hr;
for(int i=0;i<NBins-1;i++)
	{
	sprintf(targetfile,"LikelihoodDistribution/LikelihoodDistribution_%03.0lf_%03.0lf.root",PtBins[i],PtBins[i+1]);
	f=new TFile(targetfile);
	Double_t Max=0.0;
	h=(TH1F*)f->Get("Likelihood");
	hq=(TH1F*)f->Get("Likelihood_quark");
	hg=(TH1F*)f->Get("Likelihood_gluon");
	
	h->Scale(1./h->Integral("width"));
	hq->Scale(1./hq->Integral("width"));
	hg->Scale(1./hg->Integral("width"));
	
	h->SetMarkerColor(kGreen);
	hq->SetMarkerColor(kBlack);
	hg->SetMarkerColor(kRed);
		
	hr=new TH1F("deviation","deviation",100,0,1);

	for(Int_t j=0; j<=hq->GetNbinsX();j++)
		{
		Double_t Q=hq->GetBinContent(j),
			 G=hg->GetBinContent(j);
		hr->SetBinContent(j,
				1 + ((G==0)?0:G/(Q+G)) - hq->GetBinCenter(j) //1+ -- offset su 1
				);
		}
	hr->SetLineColor(kBlue);
	hr->SetMarkerColor(kBlue);
	
	Max=(Max<h->GetMaximum())?h->GetMaximum():Max;
	Max=(Max<hq->GetMaximum())?hq->GetMaximum():Max;
	Max=(Max<hg->GetMaximum())?hg->GetMaximum():Max;
	Max=(Max<hr->GetMaximum())?hr->GetMaximum():Max;
	
	TCanvas *c1=new TCanvas("Canvas","Canvas");
	c1->cd();

	hg->SetMaximum(Max*1.05);	
	hg->Draw();
	hq->Draw("same");
	h->Draw("same");
	
	hr->Draw("same");
	
	TLegend *L=new TLegend(0.35,0.70,0.55,0.85);
	L->AddEntry("Likelihood","MC mix");
	L->AddEntry("Likelihood_quark","quark");
	L->AddEntry("Likelihood_gluon","gluon");
	L->AddEntry("deviation","deviation");
	
	L->Draw();
	
	sprintf(targetfile,"LikelihoodDistribution/LikelihoodPlot_%03.0lf_%03.0lf.pdf",PtBins[i],PtBins[i+1]);
	
	//cancella la statistica
	c1->SetTitle("");
	c1->SaveAs(targetfile);
	}
}
