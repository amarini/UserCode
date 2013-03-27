

#include "TMinuit.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <cstdio>
using namespace std;

string Zselection="&& axis1_QCJet0>0 && axis2_QCJet0>0 && mZ>70 && mZ<110 && abs(ptZ-ptJet0)/(ptZ+ptJet0)<.4 && nPFCand_QC_ptCutJet>0 && ptD_QCJet0>0 && (betaStarJet0 < 0.2 * TMath::Log( nvertex - 0.67)) && (deltaPhi_jet>3.1415-0.5) ";

//----- global vars called by minuit :(
TTree *t_mc;
TTree *t_data;
string selection="ptJet0>30 && ptJet0<80 && rhoPF>0 && rhoPF<15 && .0<=abs(etaJet0) && abs(etaJet0)<2.0"+Zselection;
//string varName="QGLikelihood2012Jet0";
//string varName="QGLMLP";

string varName="QGLMLP";//QGL HISTO
string opt="CHI2 WW";
Int_t nBins=30;
Double_t xMin=0.,xMax=1.00001;
//---- internal not modify
TH1F*h_mc,*h_data;
TGraph2D *g2;

void FCN(Int_t &npar,Double_t*gin,Double_t&f,Double_t*par,Int_t flag)
{
//switch (flag)
//{
//case 1: //Init	
//	t_data->Draw( Form("%s>>h_data(%d,%lf,%lf)",varName.c_str(),nBins,xMin,xMax),selection.c_str(), "E");
//	h_data=(TH1F*)gDirectory->Get("h_data");
//	break;
//case 2:break;
//case 3:break;
//default:
	{
	//here entres parameters
	string var=Form("TMath::ATan( %f * TMath::Tan(TMath::Pi()*%s-TMath::Pi()/2.) + %f)/TMath::Pi() +0.5 ",par[0],varName.c_str(),par[1]
				);
	string sel=string("PUReWeight*eventWeight*("+selection+")");
//	cout<<"Going to Draw:"<<var<<" With Selection "<<selection<<endl;

	t_mc->Draw( Form("%s>>h_mc(%d,%lf,%lf)",var.c_str(),nBins,xMin,xMax),sel.c_str(), "E" );
	h_mc=(TH1F*)gDirectory->Get("h_mc");
	//scale
	h_mc->Scale(h_data->Integral()/h_mc->Integral());
	f=h_data->Chi2Test(h_mc,opt.c_str());
	//cout<<"CHI2="<<f<<" par[0]=" <<par[0]<<" par[1]="<<par[1]<<endl;
	}
//}
return;
}


TMinuit *gMinuit=new TMinuit(2);

int ComputeMinMinuit(){

 //Set GLobal vars
 TFile *f_data=TFile::Open("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DoubleMu-Run2012C.root");
 t_data=(TTree*)f_data->Get("tree_passedEvents");
 TFile *f_mc=TFile::Open("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_2.root");
 t_mc=(TTree*)f_mc->Get("tree_passedEvents");

gMinuit->SetFCN(FCN);

double arglist[10];
double vstrt[2];
double stp[2];
double bmin[2];
double bmax[2];

//par 4: 0.948926
//par 5: 0.097168

vstrt[0]=.9489; vstrt[1]=0.09716;
//vstrt[0]=1.0; vstrt[1]=0.;
stp[0]=0.01; stp[1]=0.01;
bmin[0]=0.9; bmax[0]=1.1;
bmin[1]=-.8; bmax[1]=.8;

//INIT 
t_data->Draw( Form("%s>>h_data(%d,%lf,%lf)",varName.c_str(),nBins,xMin,xMax),selection.c_str(), "E");
h_data=(TH1F*)gDirectory->Get("h_data");

int ierflag=0;
gMinuit->mnparm(0,"Shape",vstrt[0],stp[0],bmin[0],bmax[0],ierflag);
gMinuit->mnparm(1,"Shift",vstrt[1],stp[1],bmin[1],bmax[1],ierflag);

//SetOutput
gMinuit->SetPrintLevel(1);
//gMinuit->mnexcm("SET NOW",arglist,1,ierflag); //NO WARNING
gMinuit->SetErrorDef(1);
arglist[0]=1; //standard 1, improve minimum 2 (slower);
gMinuit->mnexcm("SET STR",arglist,1,ierflag);
//FIX PAR

arglist[0]=0.1;
gMinuit->mnexcm("SET EPS",arglist,1,ierflag);

 arglist[0] = 1;
//Double_t arglist[10];arglist[0]=0;
gMinuit->mnexcm("SCAN",arglist,1,0); //do a scan
	TGraph *gpar0=dynamic_cast<TGraph *>(gMinuit->GetPlot()->Clone("par0" ) );

 arglist[0] = 2;
gMinuit->mnexcm("SCAN",arglist,1,0); //do a scan
	TGraph *gpar1=dynamic_cast<TGraph *>(gMinuit->GetPlot()->Clone("par1") );

TCanvas *c2=new TCanvas("c2","c2",1200,600);
	c2->Divide(2);
	c2->cd(1);
	gpar0->Draw("ALP");
	c2->cd(2);
	gpar1->Draw("ALP");
TCanvas *c1=new TCanvas("c1","c1");
c1->cd();
//arglist[0]=1;gMinuit->mnexcm("FIX",arglist,1,ierflag);
//arglist[0]=500; //max iterations
//gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);

//arglist[0]=2; //standard 1, improve minimum 2 (slower);
//gMinuit->mnexcm("SET STR",arglist,1,ierflag);

//arglist[0]=500; //max iterations
//gMinuit->mnexcm("MIGRAD",arglist,1,ierflag);


//-------------------AFTER MINIMIZATION FAIL------------------
//float min0,err0;gMinuit->GetParameter(0,min0,err0);
//float min1,err1;gMinuit->GetParameter(1,min1,err1);

float min0,v0=-1;
float min1,v1=-1;

double *x,*y;
x=gpar0->GetX();y=gpar0->GetY();
for(int i=0;i<gpar0->GetN();i++){if(((v0>y[i])||(v0<0) )&&(y[i]>0)){min0=x[i];v0=y[i];}  }
x=gpar1->GetX();y=gpar1->GetY();
for(int i=0;i<gpar1->GetN();i++){if(((v1>y[i])||(v1<0) )&&(y[i]>0)){min1=x[i];v1=y[i];}  }

printf("min0=%f min1=%f, v0=%f v1=%f %d %d\n",min0,min1,v0,v1,gpar0->GetN(),gpar1->GetN());

float nstep=10;
//float stp0=.01;
//float stp1=.02;
float stp0=.03;
float stp1=.03;

//TFitResults* r0=gpar0->Fit("pol2");
//TFitResults* r1=gpar1->Fit("pol2");

TGraph2D *g=new TGraph2D(); g->SetName("Chi2"); int k=0;
for(int i=-nstep;i<=nstep;i++)
for(int j=-nstep;j<=nstep;j++)
	{
	Double_t par[4];
	par[0]=min0+i*stp0;
	par[1]=min1+j*stp1;
	Double_t f;
	FCN(4,NULL,f,par,0);
	g->SetPoint(k++,par[0],par[1],f);	
	}
TF2 *f2=new TF2("parab","[0]*(  TMath::Power((cos([1])*x+sin([1])*y+[2])/[4],2) + TMath::Power( (-sin([1])*x+cos([1])*y + [3])/[5],2)  )",0.5,1.5,0.,.5);

f2->SetParName(0,"Norm");//f2->SetParameter(0,TMath::Min(v0,v1));f2->SetParLimits(0,0,100000); //normA
f2->FixParameter(0,TMath::Min(v0,v1));
f2->SetParName(1,"alpha");f2->SetParameter(1,0); //alpha
f2->SetParName(2,"mu");f2->SetParameter(2,-min0); //mu
f2->SetParName(3,"rho");f2->SetParameter(3,-min1); //rho
f2->SetParName(4,"a");f2->SetParameter(4,.1); f2->SetParLimits(4,0,100);//a
f2->SetParName(5,"b");f2->SetParameter(5,.1); f2->SetParLimits(5,0,100);//b

g->Fit("parab");
//f2->SetParameter(0,TMath::Min(v0,v1));f2->SetParLimits(0,0,100000);
//g->Fit("parab");

double * p=f2->GetParameters();

printf("min=%.3f alpha=%.3f mu=%.3f rho=%.3f a=%.3f b=%.3f\n",p[0],p[1],p[2],p[3],p[4],p[5]);
TCanvas *c3=new TCanvas("c3","c3",800,800);
c3->Divide(2);c3->cd(1);
	//Superimposing in 3D sucks
//	double maxz=GetMaximum()*1.05;
//	double minz=g->GetMinimum()*.95;
//		g->SetMaximum(maxz);f2->SetMaximum(maxz);
//		g->SetMinimum(maxz);f2->SetMinimum(maxz);
//	double minx,maxx,miny,maxy;
	// CONT4 Changes Pad coordinates
//	TPad    *z1 = new TPad("z1","z1", 0.1, 0.1, 0.9, 0.9); z1->Draw();
//	TPad    *z2 = new TPad("z2","z2", 0.1, 0.1, 0.9, 0.9); z2->Draw();
//	z2->SetFillStyle(4000); // z2 in transparent
//	z1->cd();
g->Draw("CONT4 Z");
//	z2->cd();
//		minx=g->GetXaxis()->GetXmin();
//		maxx=g->GetXaxis()->GetXmax();
//		miny=g->GetYaxis()->GetXmin();
//		maxy=g->GetYaxis()->GetXmax();
//   	z2->Range(minx,miny,maxx,maxy);
//f2->Draw("CONT2 SAME");

c3->cd(2);
f2->Draw("CONT4");

//Found min of g
{ //<- PUOI NON IGNORARLE CINT DI MERDA?!?
float a=-1,x0=-1,y0=-1;
double * x1,*y1,*z1;
x1=g->GetX();
y1=g->GetY();
z1=g->GetZ();
for(int i=0;i<g->GetN();i++){if((z1[i]<a)||(a<0)){a=z1[i];x0=x1[i];y0=y1[i];}}
printf("min Point= %.3f %.3f %.3f\n",x0,y0,a);
}
return;
}
