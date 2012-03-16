#include "TH1F.h"
#include <stdio.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
//Fit procedure with points exclusion
#include "new/FitN.C"

inline double gammadistr_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}
inline double functionPtD_(double * x ,double*par)
{
        if((x[0]-par[0])<0)return 0;
        return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;
}
int PlotUnfoldedParams(const char *varName="nCharged",const char*Directory="Inversion",int PtMin=80,int PtMax=100)
{
TFile *f;
TCanvas *c;
TH1F*q,*g;
TH1F*rq0=new TH1F("rq0","rq0",10,0,20);
TH1F*rg0=new TH1F("rg0","rg0",10,0,20);
TH1F*rq1=new TH1F("rq1","rq1",10,0,20);
TH1F*rg1=new TH1F("rg1","rg1",10,0,20);
char fileName[1023];
TF1 *gammadistr=new TF1("gamma",gammadistr_,0,100,2);
TCanvas *c0=new TCanvas("c0","c0",800,800);
c0->Divide(3,3);
for(int i=0; i<13; i+=2)
	{
	sprintf(fileName,"../%s/%s_pt%d_%d_rho%d_%d.root",Directory,varName,PtMin,PtMax,i,i+2);
	f=TFile::Open(fileName);
	c=(TCanvas*)f->Get("c1_n3");
	q=(TH1F*)c->FindObject("Unfold_1");
	q->Scale(1./q->Integral("width"));
	c0->cd(i/2+1);
	q->Draw("P");
	TGraphErrors *Q=new TGraphErrors();int a=0;for(int k=0;k<=q->GetNbinsX();k++)if(q->GetBinError(k)>0.001){Q->SetPoint(a,q->GetBinCenter(k),q->GetBinContent(k));Q->SetPointError(a,1,q->GetBinError(k));a++;}
		gammadistr->SetParameter(1,q->GetMean());
        	gammadistr->SetParameter(0,q->GetMean()*q->GetMean()/(q->GetRMS()*q->GetRMS()));
        	Q->Fit("gamma","N Q");//N=Don't Draw
        	Q->Fit("gamma","N M Q");//N=Don't Draw M=More
		//Q->Draw("* SAME");
        	gammadistr->SetLineColor(kBlack);
        	gammadistr->SetName("gamma_quark");
        	gammadistr->DrawCopy("SAME");
        	gammadistr->SetName("gamma");
		rq0->SetBinContent(rq0->FindBin(i+1),gammadistr->GetParameter(0));
		rq0->SetBinError  (rq0->FindBin(i+1),gammadistr->GetParError(0));
		rq1->SetBinContent(rq1->FindBin(i+1),gammadistr->GetParameter(1));
		rq1->SetBinError  (rq1->FindBin(i+1),gammadistr->GetParError(1));
	g=(TH1F*)c->FindObject("Unfold_2");
	g->Scale(1./g->Integral("width"));
	c0->cd(i/2+1);
	g->Draw("P SAME");
	TGraphErrors *G=new TGraphErrors();int a=0;for(int k=0;k<=g->GetNbinsX();k++)if(g->GetBinError(k)>0.001){G->SetPoint(a,g->GetBinCenter(k),g->GetBinContent(k));G->SetPointError(a,1,g->GetBinError(k));a++;}
		gammadistr->SetParameter(1,g->GetMean());
        	gammadistr->SetParameter(0,g->GetMean()*g->GetMean()/(g->GetRMS()*g->GetRMS()));
        	G->Fit("gamma","N Q");//N=Don't Draw
        	G->Fit("gamma","N M  Q");//N=Don't Draw M=More
		//G->Draw("* SAME");
        	gammadistr->SetLineColor(kRed);
        	gammadistr->SetName("gamma_gluon");
        	gammadistr->DrawCopy("SAME");
        	gammadistr->SetName("gamma");
		rg0->SetBinContent(rg0->FindBin(i+1),gammadistr->GetParameter(0));
		rg0->SetBinError  (rg0->FindBin(i+1),gammadistr->GetParError(0));
		rg1->SetBinContent(rg1->FindBin(i+1),gammadistr->GetParameter(1));
		rg1->SetBinError  (rg1->FindBin(i+1),gammadistr->GetParError(1));
	//ora li fitto
	f->Close();
	}
rq0->SetMarkerColor(kBlack);
rq1->SetMarkerColor(kBlack);
rg0->SetMarkerColor(kRed);
rg1->SetMarkerColor(kRed);

rq0->SetMarkerStyle(20);
rq1->SetMarkerStyle(20);
rg0->SetMarkerStyle(20);
rg1->SetMarkerStyle(20);

TF1* pol1=new TF1("pol1","[0] +[1]*x",0,20);
TF1* pol0=new TF1("pol0","[0]",0,20);
TCanvas *c1=new TCanvas();
c1->Divide(2);
c1->cd(1);
rq0->Draw("P");
rg0->Draw("P SAME");
		FitN *B=new FitN();
		B->SetP(0.05);
		B->SetE(0.0001);
		B->SetN(4);
		B->SetTH1F(rq0);
		B->SetRange(0,8);
		B->SetTF1(pol0);
B->Fit();
//rq0->Fit("pol1","N Q");
pol0->SetLineColor(kBlack);pol0->DrawCopy("SAME");printf("==0 %f bq %f==\n==0 %f aq 0==\n",TMath::Log((PtMin+PtMax)/2.),pol0->GetParameter(0),TMath::Log((PtMin+PtMax)/2.));
//rg0->Fit("pol1","N Q");
B->SetTH1F(rg0);
B->Fit();
pol0->SetLineColor(kRed);pol0->DrawCopy("SAME");printf("==0 %f bg %f==\n==0 %f ag 0==\n",TMath::Log((PtMin+PtMax)/2.),pol0->GetParameter(0),TMath::Log((PtMin+PtMax)/2.));

		B->SetRange(0,12);
		B->SetTF1(pol1);
c1->cd(2);
rq1->Draw("P");
rg1->Draw("P SAME");
B->SetTH1F(rq1);
B->Fit();
//rq1->Fit("pol1","N Q");
pol1->SetLineColor(kBlack);pol1->DrawCopy("SAME");printf("==1 %f bq %f==\n==1 %f aq %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));
B->SetTH1F(rg1);
B->Fit();
//rg1->Fit("pol1","N Q");
pol1->SetLineColor(kRed);pol1->DrawCopy("SAME");printf("==1 %f bg %f==\n==1 %f ag %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));

}
