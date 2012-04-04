#include <stdio.h>

inline double gammadistr_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}

//half gamma+ offset
inline double functionPtD_(double * x ,double*par)
{
        if((x[0]-par[0])<0)return 0;
        return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;
}

//=============================================================================
int mkPlotUnfolding(const char*varName="nCharged",
int PtMin=80, int PtMax=110,int RhoMin=4,int RhoMax=6
)
//=============================================================================
{
gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetFrameFillColor(0);
//copy str
char fileName[1023];
sprintf(fileName,"../Inversion/%s_pt%d_%d_rho%d_%d.root",varName,PtMin,PtMax,RhoMin,RhoMax);

TF1 *gammadistr=new TF1("gamma",gammadistr_,0,100,2);
TF1 *functionPtD=new TF1("functionPtD",functionPtD_,0,1,3);
        functionPtD->SetParLimits(0,0.,0.4);//offset 
        functionPtD->SetParLimits(1,2,50);
        functionPtD->SetParLimits(2,0.001,0.99);

TFile *f=TFile::Open(fileName);
TCanvas *c=(TCanvas*)f->Get("c1_n3")->Clone("c1");;
TH1F*u1=(TH1F*)c->FindObject("Unfold_1");
TH1F*u2=(TH1F*)c->FindObject("Unfold_2");
TH1F*v1=(TH1F*)c->FindObject("var1");
TH1F*v2=(TH1F*)c->FindObject("var2");
TH1F*q=(TH1F*)c->FindObject("quark");
TH1F*g=(TH1F*)c->FindObject("gluon");
//PTD must be scaled with width
u1->Scale(1./u1->Integral("width"));
u2->Scale(1./u2->Integral("width"));
q->Scale(1./q->Integral("width"));
g->Scale(1./g->Integral("width"));
v1->Scale(1./v1->Integral("width"));
v2->Scale(1./v2->Integral("width"));

double Max1=u1->GetMaximum(),Max2=v1->GetMaximum();
double Min1=u1->GetMinimum(),Min2=v1->GetMinimum();
//Max1
Max1=(Max1>u2->GetMaximum())?Max1:u2->GetMaximum();Max1=(Max1>q->GetMaximum())?Max1:q->GetMaximum();Max1=(Max1>g->GetMaximum())?Max1:g->GetMaximum();
Min1=(Min1<u2->GetMinimum())?Min1:u2->GetMinimum();Min1=(Min1<q->GetMinimum())?Min1:q->GetMinimum();Min1=(Min1<g->GetMinimum())?Min1:g->GetMinimum();

//Max2
Max2=(Max2>v2->GetMaximum())?Max2:v2->GetMaximum();Max2=(Max2>q->GetMaximum())?Max2:q->GetMaximum();Max2=(Max2>g->GetMaximum())?Max2:g->GetMaximum();
Min2=(Min2<v2->GetMinimum())?Min2:v2->GetMinimum();Min2=(Min2<q->GetMinimum())?Min2:q->GetMinimum();Min2=(Min2<g->GetMinimum())?Min2:g->GetMinimum();

Max1*=1.15;
Max2*=1.15;
Min1=(Min1>0)?Min1*.85:Min1*1.15;
Min2=(Min2>0)?Min2*.85:Min2*1.15;

char xTitle[1023];
switch (varName[1])
{
case 'C':
	sprintf(xTitle,"Charged Multiplicity");
	break;
case 'N':
	sprintf(xTitle,"Neutral Multiplicity");
	break;
case 't':
	sprintf(xTitle,"PtD");
	break;
case 'G':
	sprintf(xTitle,"Q-G LD");
	q->Rebin(2);q->Scale(1./2.);
	g->Rebin(2);g->Scale(1./2.);
	u1->Rebin(2);u1->Scale(1./2.);u1->SetMarkerSize(1.6);
	u2->Rebin(2);u2->Scale(1./2.);u2->SetMarkerSize(1.6);
	v1->Rebin(2);v1->Scale(1./2.);v1->SetMarkerSize(1.6);
	v2->Rebin(2);v2->Scale(1./2.);v2->SetMarkerSize(1.6);
	break;
}

TCanvas *c1=new TCanvas("c2","c2");
u1->SetMaximum(Max1);u1->SetMinimum(Min1);
u1->GetXaxis()->SetTitle(xTitle);
u1->Draw("AXIS");
u1->Draw("AXIS X+ Y+ SAME");
q->Draw("HIST SAME");
g->Draw("HIST SAME");
u1->Draw("P SAME");
u2->Draw("P SAME");

switch (varName[1])
{
case 'C':
case 'N':
TGraphErrors *Q=new TGraphErrors();int a=0;for(int k=0;k<=u1->GetNbinsX();k++)if(u1->GetBinError(k)>0.001){Q->SetPoint(a,u1->GetBinCenter(k),u1->GetBinContent(k));Q->SetPointError(a,1,u1->GetBinError(k));a++;}
                gammadistr->SetParameter(1,u1->GetMean());
                gammadistr->SetParameter(0,u1->GetMean()*u1->GetMean()/(u1->GetRMS()*u1->GetRMS()));
                Q->Fit("gamma","N Q");//N=Don't Draw
                Q->Fit("gamma","N M Q");//N=Don't Draw M=More
                Q->Fit("gamma","N M Q");//N=Don't Draw M=More
                //Q->Draw("* SAME");
                gammadistr->SetLineColor(kBlack);
                gammadistr->SetName("gamma_quark");
                gammadistr->DrawCopy("SAME");
                gammadistr->SetName("gamma");
		printf("__QUARK__ %f %f\n",gammadistr->GetParameter(0),gammadistr->GetParameter(1));
	TGraphErrors *G=new TGraphErrors();int a=0;for(int k=0;k<=u2->GetNbinsX();k++)if(u2->GetBinError(k)>0.001){G->SetPoint(a,u2->GetBinCenter(k),u2->GetBinContent(k));G->SetPointError(a,1,u2->GetBinError(k));a++;}
                gammadistr->SetParameter(1,u2->GetMean());
                gammadistr->SetParameter(0,u2->GetMean()*u2->GetMean()/(u2->GetRMS()*u2->GetRMS()));
                G->Fit("gamma","N Q");//N=Don't Draw
                G->Fit("gamma","N M  Q");//N=Don't Draw M=More
                G->Fit("gamma","N M  Q");//N=Don't Draw M=More
                //G->Draw("* SAME");
                gammadistr->SetLineColor(kRed);
                gammadistr->SetName("gamma_gluon");
                gammadistr->DrawCopy("SAME");
                gammadistr->SetName("gamma");
		printf("__GLUON__ %f %f\n",gammadistr->GetParameter(0),gammadistr->GetParameter(1));

	break;
case 't':
	TGraphErrors *Q=new TGraphErrors();int a=0;for(int k=0;k<=u1->GetNbinsX();k++)if(u1->GetBinError(k)>0.05){Q->SetPoint(a,u1->GetBinCenter(k),u1->GetBinContent(k));Q->SetPointError(a,1./u1->GetNbinsX(),u1->GetBinError(k));a++;}
                functionPtD->SetParameter(2,u1->GetMean());
                functionPtD->SetParameter(1,(u1->GetMean())*(u1->GetMean())/(u1->GetRMS()*u1->GetRMS()));
                functionPtD->SetParameter(0,0.2);
                Q->Fit("functionPtD","N R Q");//N=Don't Draw
                Q->Fit("functionPtD","N M R Q");//N=Don't Draw M=More
                Q->Fit("functionPtD","N M R Q");//N=Don't Draw M=More
                //Q->Draw("* SAME");
                functionPtD->SetLineColor(kBlack);
                functionPtD->SetName("functionPtD_quark");
                functionPtD->DrawCopy("SAME");
                functionPtD->SetName("functionPtD");
		printf("__QUARK__ %f %f %f\n",functionPtD->GetParameter(0),functionPtD->GetParameter(1),functionPtD->GetParameter(2));
	 TGraphErrors *G=new TGraphErrors();int a=0;for(int k=0;k<=u2->GetNbinsX();k++)if(u2->GetBinError(k)>0.05){G->SetPoint(a,u2->GetBinCenter(k),u2->GetBinContent(k));G->SetPointError(a,1./u2->GetNbinsX(),u2->GetBinError(k));a++;}
                functionPtD->SetParameter(2,u2->GetMean());
                functionPtD->SetParameter(1,u2->GetMean()*u2->GetMean()/(u2->GetRMS()*u2->GetRMS()));
                functionPtD->SetParameter(0,0.2);
		if(PtMin==80)
                {
		fprintf(stderr,"\n80!!!\n");
		functionPtD->SetParameter(2,.25);
                functionPtD->SetParameter(1,10);
                functionPtD->SetParameter(0,0.05);
		}
                G->Fit("functionPtD","N Q");//N=Don't Draw
                G->Fit("functionPtD","N M Q");//N=Don't Draw M=More
                G->Fit("functionPtD","N M Q" );//N=Don't Draw M=More
                //G->Draw("* SAME");
                functionPtD->SetLineColor(kRed);
                functionPtD->SetName("functionPtD_gluon");
                functionPtD->DrawCopy("SAME");
                functionPtD->SetName("functionPtD");
		printf("__GLUON__ %f %f %f\n",functionPtD->GetParameter(0),functionPtD->GetParameter(1),functionPtD->GetParameter(2));

	break;

}

char title[1023];
char name[1023];
sprintf(title,"#splitline{%d<P_{T}[GeV]<%d}{ %d<#rho[GeV/u.a.]<%d}",PtMin,PtMax,RhoMin,RhoMax);
TLegend *L;
if(varName[0]!='Q')
	L=new TLegend(0.6,.5,0.88,.88,title);
else
	L=new TLegend(0.35,.5,0.65,.88,title);

L->SetFillColor(0);
L->AddEntry("Unfold_1","Quark (DATA)","FP");
L->AddEntry("Unfold_2","Gluon (DATA)","FP");
L->AddEntry("quark","Quark (MC)","F");
L->AddEntry("gluon","Gluon (MC)","F");
L->Draw();
sprintf(name,"../Inversion/%s_pt%d_%d_rho%d_%d_UN1.pdf",varName,PtMin,PtMax,RhoMin,RhoMax);
c1->SaveAs(name);


TCanvas *c3=new TCanvas("c3","c3");
v1->SetMaximum(Max2);v1->SetMinimum(Min2);
v1->GetXaxis()->SetTitle(xTitle);
v1->Draw("AXIS");
v1->Draw("AXIS X+ Y+ SAME");
q->Draw("HIST SAME");
g->Draw("HIST SAME");
v1->Draw("P SAME");
v2->Draw("P SAME");

if(varName[0]!='Q')
	L=new TLegend(.6,.5,0.88,.88,title);
else
	L=new TLegend(0.35,.5,0.65,.88,title);
L->SetFillColor(0);
L->AddEntry("var1","dijet","FP");
L->AddEntry("var2","#gamma jet","FP");
L->AddEntry("quark","Quark (MC)","F");
L->AddEntry("gluon","Gluon (MC)","F");
L->Draw();
sprintf(name,"../Inversion/%s_pt%d_%d_rho%d_%d_UN2.pdf",varName,PtMin,PtMax,RhoMin,RhoMax);
c3->SaveAs(name);
sprintf(name,"../Inversion/%s_pt%d_%d_rho%d_%d_UN.root",varName,PtMin,PtMax,RhoMin,RhoMax);
TFile *Out=TFile::Open(name,"RECREATE");
Out->cd();
c3->Write();
c1->Write();
Out->Close();

return 0;
}
