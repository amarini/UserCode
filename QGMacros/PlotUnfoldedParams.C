#include "TH1F.h"
#include <stdio.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
//Fit procedure with points exclusion
#include "new/FitN.C"

inline double func_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}
inline double functionPtD_(double * x ,double*par)
{
        if((x[0]-par[0])<0)return 0;
        return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;
}
int PlotUnfoldedParams(const char *varName="nCharged",const char*Directory="Inversion",int PtMin=80,int PtMax=100,const char *outTxtFileName="")
{
TFile *f;
TCanvas *c;
TH1F*rq0=new TH1F("rq0","rq0",10,0,20);
TH1F*rg0=new TH1F("rg0","rg0",10,0,20);
TH1F*rq1=new TH1F("rq1","rq1",10,0,20);
TH1F*rg1=new TH1F("rg1","rg1",10,0,20);
char fileName[1023];
TF1 *func;//new TF1("gamma",func_,0,100,2);
for(int i=0; i<20; i+=2)
	{
	sprintf(fileName,"../%s/%s_pt%d_%d_rho%d_%d_UN.root",Directory,varName,PtMin,PtMax,i,i+2);
	f=TFile::Open(fileName);
	if(f==NULL){printf("File %s does not exist. Skipped.\n",fileName);continue;}
	c=(TCanvas*)f->Get("c2");
	if(varName[0]=='n')func=c->FindObject("gamma_quark");
	else 	func=c->FindObject("functionPtD_quark");	

		rq0->SetBinContent(rq0->FindBin(i+1),func->GetParameter(0));
		rq0->SetBinError  (rq0->FindBin(i+1),func->GetParError(0));
		rq1->SetBinContent(rq1->FindBin(i+1),func->GetParameter(1));
		rq1->SetBinError  (rq1->FindBin(i+1),func->GetParError(1));

	if(varName[0]=='n') func=c->FindObject("gamma_gluon");
	else 	func=c->FindObject("functionPtD_gluon");	
		rg0->SetBinContent(rg0->FindBin(i+1),func->GetParameter(0));
		rg0->SetBinError  (rg0->FindBin(i+1),func->GetParError(0));
		rg1->SetBinContent(rg1->FindBin(i+1),func->GetParameter(1));
		rg1->SetBinError  (rg1->FindBin(i+1),func->GetParError(1));
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
FILE *fw;
if(outTxtFileName[0]!=0)fw=fopen(outTxtFileName,"a");
else fw=stdout;
//rq0->Fit("pol1","N Q");
pol0->SetLineColor(kBlack);pol0->DrawCopy("SAME");fprintf(fw,"==0 %f bq %f==\n==0 %f aq 0==\n",TMath::Log((PtMin+PtMax)/2.),pol0->GetParameter(0),TMath::Log((PtMin+PtMax)/2.));
//rg0->Fit("pol1","N Q");
B->SetTH1F(rg0);
B->Fit();
pol0->SetLineColor(kRed);pol0->DrawCopy("SAME");fprintf(fw,"==0 %f bg %f==\n==0 %f ag 0==\n",TMath::Log((PtMin+PtMax)/2.),pol0->GetParameter(0),TMath::Log((PtMin+PtMax)/2.));

		B->SetRange(0,12);
		B->SetTF1(pol1);
c1->cd(2);
rq1->Draw("P");
rg1->Draw("P SAME");
B->SetTH1F(rq1);
B->Fit();
//rq1->Fit("pol1","N Q");
pol1->SetLineColor(kBlack);pol1->DrawCopy("SAME");fprintf(fw,"==1 %f bq %f==\n==1 %f aq %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));
B->SetTH1F(rg1);
B->Fit();
//rg1->Fit("pol1","N Q");
pol1->SetLineColor(kRed);pol1->DrawCopy("SAME");fprintf(fw,"==1 %f bg %f==\n==1 %f ag %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));
return 0;
}
