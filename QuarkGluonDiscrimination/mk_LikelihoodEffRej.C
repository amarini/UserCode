#include <stdio.h>
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"

void PlotEffRej()
{
gROOT->SetBatch();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(kFALSE);
#include "PtBins.h"
char  targetfile[255];
TFile *f;
TH1F *hq,*hg;
TGraph *hr;
for(int i=0;i<NBins-1;i++)
	{
	sprintf(targetfile,"LikelihoodDistribution/LikelihoodDistribution_%03.0lf_%03.0lf.root",PtBins[i],PtBins[i+1]);
	f=new TFile(targetfile);
	
	hq=(TH1F*)f->Get("Likelihood_quark");
	hg=(TH1F*)f->Get("Likelihood_gluon");
	
	hq->Scale(1./hq->Integral("width"));
	hg->Scale(1./hg->Integral("width"));
	
	hq->SetMarkerColor(kBlack);
	hg->SetMarkerColor(kRed);
		
	hr=new TGraph(100);

	for(Int_t j=0; j<=hq->GetNbinsX();j++)
		{
		hr->SetPoint(j,
				hq->Integral(0,j,"width"),hg->Integral(0,j,"width")
				);
		}
	hr->SetLineColor(kBlue);
	hr->SetMarkerColor(kBlue);
	hr->SetLineWidth(2);	
	
	TCanvas *c1=new TCanvas("Canvas","Canvas");
	c1->cd();

	hr->Draw("ACP");
	
	sprintf(targetfile,"LikelihoodDistribution/EffRejPlot_%03.0lf_%03.0lf.pdf",PtBins[i],PtBins[i+1]);
	
	//cancella la statistica
	c1->SetTitle("");
	c1->SaveAs(targetfile);
	}
}
