#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraph2D.h"
#include "TROOT.h"
#include "TStyle.h"
/* This version wants to improve efficiency
 *
 */
using namespace std;

int PlotParams(const char *fileName="QCD_all.root",const char *variable="ptD0_quark",const char*outFileName="")
{

 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetOptStat(kFALSE);

TFile *f=TFile::Open(fileName);
TGraph2D *g=(TGraph2D*)f->Get(variable);
TCanvas *c1=new TCanvas();
g->GetXaxis()->SetTitle("P_{T}");
g->GetYaxis()->SetTitle("#rho");
g->Draw("LEGO2 A");

if(outFileName[0]!='\0')
	c1->SaveAs(outFileName);

//delete c1;
return 0;
}

