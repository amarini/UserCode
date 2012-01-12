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
/* This version wants to improve efficiency
 *
 */
using namespace std;
int getBins(double  *Bins,int nBins,double MinBin=15.0,double MaxBin=1000.,bool log=false);
int getBin(int nBins,double  *Bins,double value,double*x0=0,double*x1=0);

int PlotDiscrComp(const char *fileName="QCD_all.root",const char* fileName2="",const char *variables="ptD nCharged nNeutral")
{
double PtBins[1023];
double RhoBins[1023];
int nRhoBins=20;
int nPtBins=20;
getBins(PtBins,nPtBins,15,1000.,true);
getBins(RhoBins,nRhoBins,0,20.,false);
map< string, TH1F *> plots;
TFile *f=TFile::Open(fileName);
TFile *f2=TFile::Open(fileName2);
//general things of interest
char str[1023];
char Dir[1023];
char cut[1023];
char plotName[1023];
char VarName[1023];
const char *VariablesPointer=variables; int n;
fprintf(stderr,"Beginning var loops\n");
while(sscanf(VariablesPointer,"%s%n",VarName,&n)==1)
	{
	fprintf(stderr,"Variable=%s\n",VarName);
	VariablesPointer+=n;
	//for each bins in pt
	for(int p=0;p<nPtBins;p++)
	{
	//create the directory
	sprintf(Dir,"rhoBins_pt%.0lf_%.0lf",ceil(PtBins[p]),ceil(PtBins[p+1]));
	sprintf(str,"mkdir %s",Dir);	
	system(str);
	for(int r=0;r<nRhoBins;r++)
	{
	sprintf(plotName,"%s/%s_gluon_pt%.0lf_%.0lf_rho%.0lf",Dir,VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),ceil(RhoBins[r]));//construction of plot name
	TH1F * h0_q=(TH1F*)f->Get(plotName);
	TH1F * h1_q=(TH1F*)f2->Get(plotName);
	//quark
	sprintf(plotName,"%s/%s_quark_pt%.0lf_%.0lf_rho%.0lf",Dir,VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),ceil(RhoBins[r]));//construction of plot name
	TH1F * h0_g=(TH1F*)f->Get(plotName);
	TH1F * h1_g=(TH1F*)f2->Get(plotName);
	if((h0_q==NULL)||(h0_g==NULL)||(h1_q==NULL)||(h1_g==NULL))continue;
	
	h0_q->Scale(1./h0_q->Integral("width"));
	h0_g->Scale(1./h0_g->Integral("width"));
	h1_q->Scale(1./h1_q->Integral("width"));
	h1_g->Scale(1./h1_g->Integral("width"));
	
	float Max=0.;
	Max=(h0_q->GetMaximum()>Max)?h0_q->GetMaximum():Max;
	Max=(h1_q->GetMaximum()>Max)?h1_q->GetMaximum():Max;
	Max=(h0_g->GetMaximum()>Max)?h0_g->GetMaximum():Max;
	Max=(h1_g->GetMaximum()>Max)?h1_g->GetMaximum():Max;
	
	TCanvas *c1=new TCanvas();
	h0_q->SetLineColor(kBlack);
	h0_q->SetFillColor(kGray);
	h1_q->SetLineColor(kBlue);
	h1_q->SetMarkerColor(kBlue);
	h1_q->SetMarkerStyle(29);
	//h1_q->SetFillColor(kBlue-9);
	h0_g->SetLineColor(kRed);
	h0_g->SetMarkerColor(kRed);
	h0_g->SetMarkerStyle(29);
	h0_g->SetFillColor(kRed-9);
	h1_g->SetLineColor(kMagenta);
	h1_g->SetMarkerColor(kMagenta);
	h1_g->SetMarkerStyle(20);
	//h1_g->SetFillColor(kMagenta-9);

	h0_q->SetMaximum(Max*1.1);	
	h0_q->Draw("HIST");
	h0_g->Draw("HIST SAME");
	TH1F* tmp=(TH1F*)h0_q->Clone("tmp");tmp->SetFillColor(0);tmp->Draw("HIST SAME");
	h1_q->Draw("P SAME");
	h1_g->Draw("P SAME");

	c1->RedrawAxis();
	sprintf(str,"%s.pdf",plotName);	
	c1->SaveAs(str);
	}//loop on rho Bins
	}//loop on pt bins
	}//end variables
return 0;
}

int getBins(double  *Bins,int nBins,double MinBin,double MaxBin,bool log)
{
double incr;
if(log)
	{
	incr=TMath::Power(MaxBin/MinBin,1.0/double(nBins));
	Bins[0]=MinBin;
	Bins[nBins]=MaxBin;
	for(int i=1;i<nBins;i++)
		Bins[i]=Bins[i-1]*incr;
	}
else
	{
	incr=(MaxBin-MinBin)/nBins;
	Bins[0]=MinBin;
	Bins[nBins+1]=MaxBin;
	for(int i=1; i<nBins+1;i++)
		Bins[i]=Bins[i-1]+incr;
	}
return 0;
}
int getBin(int nBins,double  Bins[],double value,double *x0,double *x1)
{
int R=0;
//int nBins=sizeof(Bins)/sizeof(double);//?
if(value <Bins[0])return -1;
if(value >Bins[nBins])return -1;
for(R=0;R<nBins;R++)
	{
	if(Bins[R]>value)break;	
	}
R--;
if(x0) *x0=Bins[R];
if(x1) *x1=Bins[R+1];
return R;	
}
