#include "TH1F.h"
#include <stdio.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
//Fit procedure with points exclusion
#include "new/FitN.C"
#include "new/Fit2.C"
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
//TH1F*rq0=new TH1F("rq0","rq0",10,0,20);
//TH1F*rg0=new TH1F("rg0","rg0",10,0,20);
//TH1F*rq1=new TH1F("rq1","rq1",10,0,20);
//TH1F*rg1=new TH1F("rg1","rg1",10,0,20);
//TH1F*rq2=new TH1F("rq2","rq2",10,0,20);
//TH1F*rg2=new TH1F("rg2","rg2",10,0,20);
TGraphErrors*rq0=new TGraphErrors();rq0->SetName("rq0");
TGraphErrors*rg0=new TGraphErrors();rg0->SetName("rg0");
TGraphErrors*rq1=new TGraphErrors();rq1->SetName("rq1");
TGraphErrors*rg1=new TGraphErrors();rg1->SetName("rg1");
TGraphErrors*rq2=new TGraphErrors();rq2->SetName("rq2");
TGraphErrors*rg2=new TGraphErrors();rg2->SetName("rg2");
char fileName[1023];
FILE*fMeans=fopen( ("../"+string(Directory)+"/Means.txt").c_str(),"r");
printf("Opening file %s\n",("../"+string(Directory)+"/Means.txt").c_str());
if(fMeans==NULL) {printf("NO MEANS FILE\n");return 1;}
char MeansStr[1023],ch;
TF1 *func;//new TF1("gamma",func_,0,100,2);
float PtMean=-1.;
for(int i=0; i<20; i+=2)
	{
	int RhoMin=i,RhoMax=i+2;
	//Read Means.txt in order to find the correct pt bin
		printf("Opening file Means\n");
		rewind(fMeans);
		int pt0_, pt1_, rho0_,rho1_;
		PtMean=(PtMin+PtMax)/2.0,PtRMS=(PtMax-PtMin)/2.0,RhoMean=(RhoMin+RhoMax)/2.0,RhoRMS=(RhoMax-RhoMin)/2.0;	
		float ptMean=(PtMin+PtMax)/2.0,ptRMS=(PtMax-PtMin)/2.0,rhoMean=(RhoMin+RhoMax)/2.0,rhoRMS=(RhoMax-RhoMin)/2.0;	
		while (fscanf(fMeans,"%d:%d %d:%d %f %f %f %f",&pt0_,&pt1_,&rho0_,&rho1_,&ptMean,&ptRMS,&rhoMean,&rhoRMS)!=EOF){
			if((pt0_==PtMin)&&(pt1_==PtMax)&&(rho0_==RhoMin)&&(rho1_==RhoMax)){
					PtMean=ptMean;
					PtRMS=ptRMS;
					RhoMean=rhoMean;
					RhoRMS=rhoRMS;
				};
			};
		printf("DONE\n");
	//
	sprintf(fileName,"../%s/%s_pt%d_%d_rho%d_%d_UN.root",Directory,varName,PtMin,PtMax,i,i+2);
	f=TFile::Open(fileName);
	if(f==NULL){printf("File %s does not exist. Skipped.\n",fileName);continue;}
	c=(TCanvas*)f->Get("c2");
	if(varName[0]=='n')func=c->FindObject("gamma_quark");
	else 	func=c->FindObject("functionPtD_quark");	

		rq0->SetPoint(rq0->GetN(),RhoMean,func->GetParameter(0));
		rq0->SetPointError  (rq0->GetN()-1,RhoRMS,func->GetParError(0));
		rq1->SetPoint(rq1->GetN(),RhoMean,func->GetParameter(1));
		rq1->SetPointError  (rq1->GetN()-1,RhoRMS,func->GetParError(1));
	if(varName[0]!='n'){
		rq2->SetPoint(rq2->GetN(),RhoMean,func->GetParameter(2));
		rq2->SetPointError  (rq2->GetN()-1,RhoRMS,func->GetParError(2));
		}//

	if(varName[0]=='n') func=c->FindObject("gamma_gluon");
	else 	func=c->FindObject("functionPtD_gluon");	

		rg0->SetPoint(rg0->GetN(),RhoMean,func->GetParameter(0));
		rg0->SetPointError  (rg0->GetN()-1,RhoRMS,func->GetParError(0));
		rg1->SetPoint(rg1->GetN(),RhoMean,func->GetParameter(1));
		rg1->SetPointError  (rg1->GetN()-1,RhoRMS,func->GetParError(1));
	if(varName[0]!='n'){
		rg2->SetPoint(rg2->GetN(),RhoMean,func->GetParameter(2));
		rg2->SetPointError  (rg2->GetN()-1,RhoRMS,func->GetParError(2));
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
rq0->Draw("A P");
rg0->Draw("P SAME");

Fit2 *D=new Fit2();
D->SetTGRAPH(rq0,0);
D->SetTGRAPH(rg0,1);
//D->SetTH1F(rq0,0);
//D->SetTH1F(rg0,1);
std::vector<float> *v=D->Fit();
pol1->SetParameter(0,v->at(0));
pol1->SetParameter(1,v->at(1));

FILE *fw;
if(outTxtFileName[0]!=0)fw=fopen(outTxtFileName,"a");
else fw=stdout;

pol1->SetLineColor(kBlack);pol1->DrawCopy("SAME");fprintf(fw,"==%s:0 %f bq %f==\n==%s:0 %f aq %f==\n",varName,TMath::Log(PtMean),pol1->GetParameter(0),varName,TMath::Log(PtMean),pol1->GetParameter(1));

pol1->SetParameter(0,v->at(2));
pol1->SetParameter(1,v->at(1));
delete v;
pol1->SetLineColor(kRed);pol1->DrawCopy("SAME");fprintf(fw,"==%s:0 %f bg %f==\n==%s:0 %f ag %f==\n",varName,TMath::Log(PtMean),pol1->GetParameter(0),varName,TMath::Log(PtMean),pol1->GetParameter(1));

c1->cd(2);
D->SetTGRAPH(rq1,0);
D->SetTGRAPH(rg1,1);
v=D->Fit();
pol1->SetParameter(0,v->at(0));
pol1->SetParameter(1,v->at(1));

rq1->Draw("A P");
rg1->Draw("P SAME");

pol1->SetLineColor(kBlack);pol1->DrawCopy("SAME");fprintf(fw,"==%s:1 %f bq %f==\n==%s:1 %f aq %f==\n",varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));

pol1->SetParameter(0,v->at(2));
pol1->SetParameter(1,v->at(1));
delete v;


pol1->SetLineColor(kRed);pol1->DrawCopy("SAME");fprintf(fw,"==%s:1 %f bg %f==\n==%s:1 %f ag %f==\n",varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));

if(varName[0]!='n'){
c1->cd(3);
rq2->Draw("A P");
rg2->Draw("P SAME");
D->SetTGRAPH(rq1,0);
D->SetTGRAPH(rg1,1);
v=D->Fit();

pol1->SetParameter(0,v->at(0));
pol1->SetParameter(1,v->at(1));
pol1->SetLineColor(kBlack);pol1->DrawCopy("SAME");fprintf(fw,"==%s: 1 %f bq %f==\n==%s:1 %f aq %f==\n",varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));

pol1->SetParameter(0,v->at(2));
pol1->SetParameter(1,v->at(1));

pol1->SetLineColor(kRed);pol1->DrawCopy("SAME");fprintf(fw,"==%s:1 %f bg %f==\n==%s:1 %f ag %f==\n",varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(0),varName,TMath::Log((PtMin+PtMax)/2.),pol1->GetParameter(1));
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
