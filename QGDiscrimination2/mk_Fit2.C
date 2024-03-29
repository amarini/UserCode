#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "PtBins.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"

#include "QGLikelihoodCalculator.h"
#include "QGLikelihoodCalculator.C"

//inline double gammadistr_(double* x, double* par)
//{
//	return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;		
//}
//
//inline double functionPtD_(double * x ,double*par)
//{
//	return TMath::Exp ( (x[0]-par[0])/par[1] *(
//						-(x[0]-par[0])/par[1]+ TMath::Sqrt( TMath::Power( (x[0]-par[0])/par[1],2) +par[4] ) -par[2]
//						) 
//			
//	) * par[3]
//		* par[1]*(TMath::Sqrt(2*TMath::Pi())/2 +1./par[2] ); //normalizzazione a meta' (quasi 1 a parte per par3 che e' un'integrazione di una gauss)
//}


void mk_Fit2()
{
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetFrameFillColor(0);


char filename[1023];
char canvasname[1023];
char histoname[1023];

FILE *R=fopen("/tmp/amarini/fitresults.txt","w");
TString VarNames[] ={"ncharged","nneutral","PtD"};//,"rRMS"};
//double (*gammadistr_1)(double *, double *)= &(QGLikelihoodCalculator::gammadistr_);
//double (*functionPtD_1)(double *, double *)= &(QGLikelihoodCalculator::functionPtD_);
TF1 *gammadistr=new TF1("gamma",QGLikelihoodCalculator::gammadistr_,0,100,2);
TF1 *functionPtD=new TF1("functionPtD",QGLikelihoodCalculator::functionPtD_,0,1,5);
TFile *F=TFile::Open("/home/Giochi/flat_Z2_nopileup_new.root");
TTree *t=(TTree*)F->Get("demo/t");

for(int j=0; j<3;j++)
{
//il set dei parametri e' messo qui in modo da 'seguire' i risultati dei fit precedenti
gammadistr->SetParLimits(0,1,20);
gammadistr->SetParLimits(1,1,50);

int count=0;
TGraph*Mean_quark=new TGraph();//alpha*alpha = par[0]
TGraph*Mean_gluon=new TGraph();
TGraph*alpha_quark=new TGraph();//alpha= par[1]/par[0]
TGraph*alpha_gluon=new TGraph();
TGraph *PtD_quark[5];
TGraph *PtD_gluon[5];
char str[1023];
for(int k=0; k<5;k++)
	{
	PtD_quark[k]=new TGraph();
	PtD_gluon[k]=new TGraph();
	sprintf(str,"PtD%d_quark",k);
	PtD_quark[k]->SetName(str);
	sprintf(str,"PtD%d_gluon",k);
	PtD_gluon[k]->SetName(str);
	}

for(int i=0; i<NBins-1;i++)
{
//sprintf(filename,"nopileup/VarDistribution/VarDistribution_%03.0lf_%03.0lf.root",PtBins[i],PtBins[i+1]);
//fprintf(stderr,"Opening file: %s\n",filename);
//TFile *f=TFile::Open(filename);
char cut[1023];
char sel[1023];

sprintf(histoname,"%s_quark",VarNames[j].Data());

if( (VarNames[j].Data())[0] =='P') sprintf(sel,"%s>>%s(100,0,1)",VarNames[j].Data(),histoname);
else sprintf(sel,"%s>>%s(100,0,100)",VarNames[j].Data(),histoname);

sprintf(cut,"%lf<jtpt&&jtpt<%lf &&abs(jteta)<2.0 && abs(pdgid)<5",PtBins[i],PtBins[i+1]);
t->Draw(sel,cut);
//fprintf(stderr,"Getting Histo:%s\n",histoname);
TH1F *Histo_quark=(TH1F*)gDirectory->Get(histoname)->Clone();
sprintf(histoname,"%s_gluon",VarNames[j].Data());

if( (VarNames[j].Data())[0] =='P') sprintf(sel,"%s>>%s(100,0,1)",VarNames[j].Data(),histoname);
else sprintf(sel,"%s>>%s(100,0,100)",VarNames[j].Data(),histoname);
sprintf(cut,"%lf<jtpt&&jtpt<%lf &&abs(jteta)<2.0 && abs(pdgid)==21",PtBins[i],PtBins[i+1]);
t->Draw(sel,cut);
TH1F *Histo_gluon=(TH1F*)gDirectory->Get(histoname)->Clone();

sprintf(sel,"jtpt>>pthisto(50,%lf,%lf)",PtBins[i],PtBins[i+1]);
sprintf(cut,"%lf<jtpt&&jtpt<%lf &&abs(jteta)<2.0 ",PtBins[i],PtBins[i+1]);
t->Draw(sel,cut);
TH1F *Histo_pt=(TH1F*)gDirectory->Get("pthisto")->Clone();
double ptmean=Histo_pt->GetMean();


Histo_quark->Scale(1./Histo_quark->Integral("width"));
Histo_gluon->Scale(1./Histo_gluon->Integral("width"));


gammadistr->SetParameter(0,Histo_quark->GetRMS()*Histo_quark->GetRMS()/Histo_quark->GetMean());
gammadistr->SetParameter(1,Histo_quark->GetMean());
printf("Mean %lf - alpha %lf \n ",Histo_quark->GetMean(),Histo_quark->GetRMS()*Histo_quark->GetRMS()/Histo_quark->GetMean());

if( (VarNames[j].Data())[0] =='r')
	{
	gammadistr->SetParameter(0,6);//shape
	gammadistr->SetParameter(1,0.005);//media
	}
	
	functionPtD->SetParameter(0,0.4);
	functionPtD->SetParLimits(0,0.1,0.9);//media

	functionPtD->SetParameter(1,0.5);//sigma gaus
	functionPtD->SetParLimits(1,0.,1.);

	functionPtD->SetParameter(2,5);
	functionPtD->SetParLimits(2,0.,20);

	functionPtD->SetParameter(3,0.9);
	functionPtD->SetParLimits(3,0.,500);//fattore di normalizzazione ~1

	functionPtD->SetParameter(4,2);
	functionPtD->SetParLimits(4,0.1,50.);


if( (VarNames[j].Data())[0] =='P')
{
	Histo_quark->Fit("functionPtD","N");
	Histo_quark->Fit("functionPtD","NM");
}
else
	{
	Histo_quark->Fit("gamma","N");//N=Don't Draw
	Histo_quark->Fit("gamma","NM");//N=Don't Draw
	}

fprintf(R,"%s_quark %03.0lf_%03.0lf: %lf %lf %lf %lf\n",VarNames[j].Data(),
							PtBins[i],
							PtBins[i+1],
							gammadistr->GetParameter(0),
							gammadistr->GetParameter(1),
							gammadistr->GetParError(0),
							gammadistr->GetParError(1));
Mean_quark->SetPoint(count,ptmean,gammadistr->GetParameter(1));
alpha_quark->SetPoint(count,ptmean,gammadistr->GetParameter(0));

for(int k=0;k<5;k++)PtD_quark[k]->SetPoint(count,ptmean,functionPtD->GetParameter(k));

TCanvas *c1 =new TCanvas();
//c1->SetLogy();
Histo_quark->Draw();
gammadistr->SetLineColor(kBlack);
functionPtD->SetLineColor(kBlack);

if( (VarNames[j].Data())[0] =='P')
	functionPtD->DrawCopy("SAME C");
else
	gammadistr->DrawCopy("SAME C");

if( (VarNames[j].Data())[0] =='P')
	{
	Histo_gluon->Fit("functionPtD","N");
	Histo_gluon->Fit("functionPtD","N");
	Histo_gluon->Fit("functionPtD","NM");
	}
else
	{
	Histo_gluon->Fit("gamma","N");
	Histo_gluon->Fit("gamma","NM");
	}
	
fprintf(R,"%s_gluon %03.0lf_%03.0lf: %lf %lf %lf %lf\n",VarNames[j].Data(),
							PtBins[i],
							PtBins[i+1],
							gammadistr->GetParameter(0),
							gammadistr->GetParameter(1),
							gammadistr->GetParError(0),
							gammadistr->GetParError(1));
Mean_gluon->SetPoint(count,ptmean,gammadistr->GetParameter(1));
alpha_gluon->SetPoint(count,ptmean,gammadistr->GetParameter(0));
for(int k=0;k<5;k++)PtD_gluon[k]->SetPoint(count,ptmean,functionPtD->GetParameter(k));
count++;

Histo_gluon->SetLineColor(kRed);
Histo_gluon->Draw("SAME");
gammadistr->SetLineColor(kRed);
functionPtD->SetLineColor(kRed);
if( (VarNames[j].Data())[0] =='P')
	functionPtD->DrawCopy("SAME C");
else
	gammadistr->DrawCopy("SAME C");

//c1->SetLogy();
sprintf(canvasname,"/tmp/amarini/%s_%03.0lf_%03.0lf.png",VarNames[j].Data(),PtBins[i],PtBins[i+1]);
c1->SaveAs(canvasname);
delete c1;
//f->Close();
}



TCanvas *c1 =new TCanvas();
c1->SetLogx();


Mean_quark->SetName("mean_quark");
Mean_gluon->SetName("mean_gluon");
alpha_quark->SetName("alpha_quark");
alpha_gluon->SetName("alpha_gluon");

Mean_quark->SetLineColor(kBlack);
Mean_quark->SetMarkerColor(kBlack);
Mean_quark->SetMarkerStyle(20);
Mean_gluon->SetLineColor(kRed);
Mean_gluon->SetMarkerColor(kRed);
Mean_gluon->SetMarkerStyle(20);

alpha_quark->SetLineColor(kBlue);
alpha_quark->SetMarkerColor(kBlue);
alpha_quark->SetMarkerStyle(20);
alpha_gluon->SetLineColor(kOrange);
alpha_gluon->SetMarkerColor(kOrange);
alpha_gluon->SetMarkerStyle(20);

Mean_gluon->SetMaximum(25);
Mean_gluon->TGraph::SetMinimum(0.0);
Mean_gluon->GetXaxis()->SetMoreLogLabels();
Mean_gluon->GetXaxis()->SetNoExponent();
Mean_gluon->GetXaxis()->SetTitle("P_{T} [GeV/c]");
Mean_gluon->GetYaxis()->SetTitle("Fit Parameters");

Mean_gluon->Draw("ALP");
Mean_quark->Draw("LP SAME");
alpha_quark->Draw("LP SAME");
alpha_gluon->Draw("LP SAME");

//sprintf(str,"%s",VarNames[j].Data());
TLegend *L=new TLegend(0.2,0.6,0.45,0.85,VarNames[j].Data());
L->SetFillColor(0);
L->AddEntry("mean_quark","#beta quark","P");
L->AddEntry("mean_gluon","#beta gluon","P");
L->AddEntry("alpha_quark","#alpha quark","P");
L->AddEntry("alpha_gluon","#alpha gluon","P");
L->Draw();

sprintf(canvasname,"/tmp/amarini/%s_fitresults.png",VarNames[j].Data());
c1->SaveAs(canvasname);
sprintf(canvasname,"/tmp/amarini/%s_fitresults.root",VarNames[j].Data());
TFile *fr=new TFile(canvasname,"RECREATE");
if( (VarNames[j].Data())[0] =='P')
	{
		fr->cd();
		for(int k=0;k<5;k++)
		{
			PtD_quark[k]->Write();
			PtD_gluon[k]->Write();
			functionPtD->Write();
		}
	}
else
	{
		fr->cd();
		alpha_quark->Write();
                alpha_gluon->Write();
                Mean_quark->Write();
                Mean_gluon->Write();

		gammadistr->Write();
	}
fr->Close();
delete alpha_quark;
delete alpha_gluon;
delete Mean_quark;
delete Mean_gluon;
delete c1;

}
return ;
}
