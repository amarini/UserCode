#include "TH1F.h"
#include <stdio.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
//Fit procedure with points exclusion
#include "new/FitN.C"
#include "new/GetChar.C" //get a char without changing position

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
TH1F*rq2=new TH1F("rq2","rq2",10,0,20);
TH1F*rg2=new TH1F("rg2","rg2",10,0,20);
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
	if(varName[0]!='n'){
		rq2->SetBinContent(rq2->FindBin(i+1),func->GetParameter(2));
		rq2->SetBinError  (rq2->FindBin(i+1),func->GetParError(2));
		}//

	if(varName[0]=='n') func=c->FindObject("gamma_gluon");
	else 	func=c->FindObject("functionPtD_gluon");	

		rg0->SetBinContent(rg0->FindBin(i+1),func->GetParameter(0));
		rg0->SetBinError  (rg0->FindBin(i+1),func->GetParError(0));
		rg1->SetBinContent(rg1->FindBin(i+1),func->GetParameter(1));
		rg1->SetBinError  (rg1->FindBin(i+1),func->GetParError(1));
	if(varName[0]!='n'){
		rg2->SetBinContent(rg2->FindBin(i+1),func->GetParameter(2));
		rg2->SetBinError  (rg2->FindBin(i+1),func->GetParError(2));
		}//
	//ora li fitto
	f->Close();
	}
rq0->SetMarkerColor(kBlack);
rq1->SetMarkerColor(kBlack);
rq2->SetMarkerColor(kBlack);
rg0->SetMarkerColor(kRed);
rg1->SetMarkerColor(kRed);
rg2->SetMarkerColor(kRed);

rq0->SetMarkerStyle(20);
rq1->SetMarkerStyle(20);
rq2->SetMarkerStyle(20);
rg0->SetMarkerStyle(20);
rg1->SetMarkerStyle(20);
rg2->SetMarkerStyle(20);

TF1* pol1=new TF1("pol1","[0] +[1]*x",0,20);
TF1* pol0=new TF1("pol0","[0]",0,20);
TCanvas *c1=new TCanvas();
if(varName[0]=='n')
	c1->Divide(2);
else c1->Divide(3);
c1->cd(1);
rq0->Draw("P");
rg0->Draw("P SAME");
		FitN *B=new FitN();
		B->SetP(0.05); //default 0.05
		B->SetE(0.0001);
		B->SetN(5);
		B->SetTH1F(rq0);
		B->SetRange(0,20);
		B->SetTF1(pol1);
B->Fit();
FILE *fw;
if(outTxtFileName[0]!=0)fw=fopen(outTxtFileName,"a");
else fw=stdout;
//rq0->Fit("pol1","N Q");
pol1->SetLineColor(kBlack);pol1->DrawCopy("SAME");fprintf(fw,"==0 %f bq %f==\n==0 %f aq %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));

//rg0->Fit("pol1","N Q");
B->SetTH1F(rg0);
B->Fit();
pol1->SetLineColor(kRed);pol1->DrawCopy("SAME");fprintf(fw,"==0 %f bg %f==\n==0 %f ag %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));

		B->SetRange(1,18);
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

if(varName[0]!='n'){
c1->cd(3);
rq2->Draw("P");
rg2->Draw("P SAME");
B->SetTH1F(rq2);
B->Fit();
//rq2->Fit("pol1","N Q");
pol1->SetLineColor(kBlack);pol1->DrawCopy("SAME");fprintf(fw,"==1 %f bq %f==\n==1 %f aq %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));
B->SetTH1F(rg2);
B->Fit();
//rg2->Fit("pol1","N Q");
pol1->SetLineColor(kRed);pol1->DrawCopy("SAME");fprintf(fw,"==1 %f bg %f==\n==1 %f ag %f==\n",TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));
}//end 3rd params plot

//Getting old parameters
FILE *fr=fopen( (string("../Fit/")+varName+".txt").c_str(),"r");
if(fr==NULL)return 0;
TGraph *Q_Fit0=new TGraph();Q_Fit0->SetName("Q_Fit0");
TGraph *G_Fit0=new TGraph();G_Fit0->SetName("G_Fit0");
TGraph *Q_Fit1=new TGraph();Q_Fit1->SetName("Q_Fit1");
TGraph *G_Fit1=new TGraph();G_Fit1->SetName("G_Fit1");
TGraph *Q_Fit2=new TGraph();Q_Fit2->SetName("Q_Fit2");
TGraph *G_Fit2=new TGraph();G_Fit2->SetName("G_Fit2");

//skip lines that begin with [ or { -> useless
char ch;
int count=0;
while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
        {
        ch='\0';while(ch!='\n'){fscanf(fr,"%c",&ch);fprintf(stderr,"%c",ch);}
        } //[gluon] {bla}
while((GetChar(fr) != '['))
	{
	float pt0,pt1,rho0,rho1,xmin,xmax; int n;	
	float par[10];
	fscanf(fr,"%f %f %f %f %d %f %f",&pt0,&pt1,&rho0,&rho1,&n,&xmin,&xmax);
	for(int i=0;i<n-2;i++) fscanf(fr,"%f",&par[i]);
	
	if( (PtMin<=pt0 && pt0<=PtMax) || (PtMin<=pt1&&pt1<=PtMax)) //if all the bin intervals in mc is inside
		{
		Q_Fit0->SetPoint(count, (rho1+rho0)/2.0,par[0]);
		Q_Fit1->SetPoint(count, (rho1+rho0)/2.0,par[1]);
		if(varName[0]!='n')Q_Fit2->SetPoint(count, (rho1+rho0)/2.0,par[2]);
		count++;
		}
	}
printf("Points=%d\n",count);
count=0;
while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
        {
        ch='\0';while(ch!='\n'){fscanf(fr,"%c",&ch);fprintf(stderr,"%c",ch);}
        } //[gluon] {bla}
while((GetChar(fr) != '[') && (GetChar(fr)!=EOF))
	{
	float pt0,pt1,rho0,rho1,xmin,xmax; int n;	
	float par[10];
	fscanf(fr,"%f %f %f %f %d %f %f",&pt0,&pt1,&rho0,&rho1,&n,&xmin,&xmax);
	for(int i=0;i<n-2;i++) fscanf(fr,"%f",&par[i]);
	
	if( (PtMin<=pt0 && pt0<=PtMax) || (PtMin<=pt1&&pt1<=PtMax)) //if all the bin intervals in mc is inside
		{
		G_Fit0->SetPoint(count,(rho1+rho0)/2.0,par[0]);
		G_Fit1->SetPoint(count,(rho1+rho0)/2.0,par[1]);
		if(varName[0]!='n')G_Fit2->SetPoint(count,(rho1+rho0)/2.0,par[2]);
		count++;
		}
	}
c1->cd(1);
Q_Fit0->SetMarkerColor(kGray+2);
Q_Fit0->SetMarkerStyle(29);
G_Fit0->SetMarkerColor(kRed+2);
G_Fit0->SetMarkerStyle(29);
Q_Fit0->Draw("P SAME");
G_Fit0->Draw("P SAME");
c1->cd(2);
Q_Fit1->SetMarkerColor(kGray+2);
Q_Fit1->SetMarkerStyle(29);
G_Fit1->SetMarkerColor(kRed+2);
G_Fit1->SetMarkerStyle(29);
Q_Fit1->Draw("P SAME");
G_Fit1->Draw("P SAME");
if(varName[0]!='n'){
c1->cd(3);
Q_Fit2->SetMarkerColor(kGray+2);
Q_Fit2->SetMarkerStyle(29);
G_Fit2->SetMarkerColor(kRed+2);
G_Fit2->SetMarkerStyle(29);
Q_Fit2->Draw("P SAME");
G_Fit2->Draw("P SAME");
}
return 0;
}
