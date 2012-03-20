#include "new/ComputeMixture.C"
#include "new/PtBins.h"
#include "new/PtBins.C"
int PlotMixture()
{

TGraphErrors *Ph_q=new TGraphErrors();Ph_q->SetName("ph_q");
TGraphErrors *Di_q=new TGraphErrors();Di_q->SetName("di_q");
Di_q->SetMarkerColor(kBlue+2);
Di_q->SetMarkerStyle(30);
Di_q->SetMarkerSize(1.2);
Ph_q->SetMarkerColor(kBlue+2);
Ph_q->SetMarkerStyle(30);
Ph_q->SetMarkerSize(1.2);

TGraphErrors *Ph_g=new TGraphErrors();Ph_g->SetName("ph_g");
TGraphErrors *Di_g=new TGraphErrors();Di_g->SetName("di_g");
Di_g->SetMarkerColor(kRed+2);
Di_g->SetMarkerStyle(29);
Di_g->SetMarkerSize(1.2);
Ph_g->SetMarkerColor(kRed+2);
Ph_g->SetMarkerStyle(29);
Ph_g->SetMarkerSize(1.2);

Ph_q->SetFillColor(0);
Di_q->SetFillColor(0);
Ph_g->SetFillColor(0);
Di_g->SetFillColor(0);
int count=0;
//Getting PtBinning
double PtBins[25];
double RhoBins[25];
Bins *A=new Bins();
A->getBins_int(18,PtBins,20,1000,true);
PtBins[18]=3500;
A->getBins_int(21,RhoBins,0,20,false);
//
for(int i=0;i<18;i++)
	{
	double x,ex,y,ey;
	ComputeMixture("~/Documents/QGDiscriminator/QG/Omog_DiJet_QCD_*.root",PtBins[i],PtBins[i+1],0,20,"omog",&x,&ex,&y,&ey);
	Di_q->SetPoint(count,(PtBins[i]+PtBins[i+1]) /2.0,x);
	Di_q->SetPointError(count,(PtBins[i+1]-PtBins[i]) /2.0,ex);
	Di_g->SetPoint(count,(PtBins[i]+PtBins[i+1]) /2.0,y);
	Di_g->SetPointError(count,(PtBins[i+1]-PtBins[i]) /2.0,ey);
	ComputeMixture("~/Documents/QGDiscriminator/QG/Omog_QGStudies_*_Summer11.root",PtBins[i],PtBins[i+1],0,20,"omog",&x,&ex,&y,&ey);
	Ph_q->SetPoint(count,(PtBins[i]+PtBins[i+1]) /2.0,x);
	Ph_q->SetPointError(count,(PtBins[i+1]-PtBins[i]) /2.0,ex);
	Ph_g->SetPoint(count,(PtBins[i]+PtBins[i+1]) /2.0,y);
	Ph_g->SetPointError(count,(PtBins[i+1]-PtBins[i]) /2.0,ey);
	count++;
	}

TCanvas *c1=new TCanvas("c1","c1",800,800);
c1->SetLogx();
Di_q->SetMinimum(0);
Di_q->SetMaximum(1);
Di_q->Draw("AP");
Di_g->Draw("P SAME");
TLegend *L=new TLegend(0.8,.12,.89,.25);
L->AddEntry("di_q","Dijet Quark","PF");
L->AddEntry("di_g","Dijet Quark","PF");
L->Draw();

TCanvas *c2=new TCanvas("c2","c2",800,800);
c2->SetLogx();
Ph_q->SetMinimum(0);
Ph_q->SetMaximum(1);
Ph_q->Draw("AP");
Ph_g->Draw("P SAME");
TLegend *L=new TLegend(0.8,.12,.89,.25);
L->AddEntry("ph_q","#gamma Jet Quark","PF");
L->AddEntry("ph_g","#gamma Jet Gluon","PF");
L->Draw();
return 0;
}
