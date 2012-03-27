#include "new/ComputeMixture.C"
#include "new/PtBins.h"
#include "new/PtBins.C"
#include "new/CMSLatex.h"
#include "new/CMSLatex.C"
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

TGraphErrors *EM_q=new TGraphErrors();EM_q->SetName("EM_q"); //EM enriched
TGraphErrors *EM_g=new TGraphErrors();EM_g->SetName("EM_g");
EM_q->SetMarkerColor(kOrange+2);
EM_q->SetMarkerStyle(24);
EM_q->SetMarkerSize(1.2);
EM_g->SetMarkerColor(kGreen+2);
EM_g->SetMarkerStyle(24);
EM_g->SetMarkerSize(1.2);

Ph_q->SetFillColor(0);
Di_q->SetFillColor(0);
EM_q->SetFillColor(0);
Ph_g->SetFillColor(0);
Di_g->SetFillColor(0);
EM_g->SetFillColor(0);


TGraphErrors *Di_up=new TGraphErrors();Di_up->SetName("di_up");
TGraphErrors *Di_dn=new TGraphErrors();Di_dn->SetName("di_dn");
Di_up->SetLineColor(kBlue+2);
Di_dn->SetLineColor(kBlue+2);
Di_up->SetLineStyle(2);
Di_dn->SetLineStyle(3);
Di_up->SetLineWidth(2);
Di_dn->SetLineWidth(2);
Di_up->SetFillColor(0);
Di_dn->SetFillColor(0);

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
	//Dijet
	ComputeMixture("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_QCD_*.root",PtBins[i],PtBins[i+1],0,30,"omog",&x,&ex,&y,&ey);
	Di_q->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,x);
	Di_q->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ex);
	Di_g->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,y);
	Di_g->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ey);
	
	//Photon
	ComputeMixture("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_*_Summer11.root",PtBins[i],PtBins[i+1],0,30,"omog",&x,&ex,&y,&ey);
	Ph_q->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,x);
	Ph_q->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ex);
	Ph_g->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,y);
	Ph_g->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ey);
	
	//EM	
	ComputeMixture("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_QGStudies_*EM*_Summer11.root",PtBins[i],PtBins[i+1],0,30,"omog",&x,&ex,&y,&ey);
	EM_q->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,x);
	EM_q->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ex);
	EM_g->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,y);
	EM_g->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ey);
	//Di UP
	ComputeMixture("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_QCD_*.root",PtBins[i],PtBins[i+1],10,30,"omog",&x,&ex,&y,&ey);
	Di_up->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,x);
	Di_up->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ex);
	//Di Down
	ComputeMixture("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/Omog_DiJet_QCD_*.root",PtBins[i],PtBins[i+1],0,5,"omog",&x,&ex,&y,&ey);
	Di_dn->SetPoint(i,(PtBins[i]+PtBins[i+1]) /2.0,x);
	Di_dn->SetPointError(i,(PtBins[i+1]-PtBins[i]) /2.0,ex);
	
	
	}

double a,b;

Di_dn->GetPoint(0,a,b);
while(!(-.1<b && b<1.1)){Di_dn->RemovePoint(0);Di_dn->GetPoint(0,a,b);}
Di_up->GetPoint(0,a,b);
while(!(-.1<b && b<1.1)){Di_up->RemovePoint(0);Di_up->GetPoint(0,a,b);}


TCanvas *c1=new TCanvas("c1","c1",800,800);
c1->SetLogx();
Di_q->SetMinimum(0);
Di_q->SetMaximum(1);
Di_q->GetXaxis()->SetTitle("P_{T}");
Di_q->GetYaxis()->SetTitle("Fraction");
Di_q->Draw("AP");
Di_g->Draw("P SAME");
EM_q->Draw("P SAME");
EM_g->Draw("P SAME");
Di_up->Draw("L SAME");
Di_dn->Draw("L SAME");
TLegend *L=new TLegend(0.3,.12,.7,.3);
L->AddEntry("di_q","Dijet Quark","PF");
L->AddEntry("di_g","Dijet Gluon","PF");
L->AddEntry("EM_q","EM Quark","PF");
L->AddEntry("EM_g","EM Gluon","PF");
L->AddEntry("di_up","Dijet Quark (#rho>10)","F");
L->AddEntry("di_dn","Dijet Quark (#rho<5)","F");
L->Draw();
CMSLatex::DrawSimulation();
c1->SaveAs("../Inversion/Mixture1.root");

TCanvas *c2=new TCanvas("c2","c2",800,800);
c2->SetLogx();
Ph_q->SetMinimum(0);
Ph_q->SetMaximum(1);
Ph_q->GetXaxis()->SetTitle("P_{T}");
Ph_q->GetYaxis()->SetTitle("Fraction");
Ph_q->Draw("AP");
Ph_g->Draw("P SAME");
TLegend *L=new TLegend(0.3,.35,.7,.65);
L->AddEntry("ph_q","#gamma Jet Quark","PF");
L->AddEntry("ph_g","#gamma Jet Gluon","PF");
L->Draw();
CMSLatex::DrawSimulation();
c2->SaveAs("../Inversion/Mixture2.root");
return 0;
}
