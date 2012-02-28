#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"


#include "stdio.h"
#include "stdlib.h"

#include "PtBins.h"

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


void mk_Fit2(const char * fileName,const char*outFileName="",const char * txtFileName="tmp.txt")
{
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetFrameFillColor(0);


char cmd[1023];

TString VarNames[] ={"nCharged","nNeutral","ptD"};//,"rRMS"};
TF1 *gammadistr=new TF1("gamma",gammadistr_,0,100,2);
TF1 *functionPtD=new TF1("functionPtD",functionPtD_,0,1,3);//

fprintf(stderr,"Opening Files\n");
TFile *f=TFile::Open(fileName);

TFile *F=TFile::Open(outFileName,"RECREATE");

FILE *fw=fopen(txtFileName,"w");


fprintf(stderr,"Creating Graphs\n");
F->cd();
TGraph2D* nCharged_quark[2];
TGraph2D* nCharged_gluon[2];
TGraph2D* nNeutral_quark[2];
TGraph2D* nNeutral_gluon[2];
TGraph2D* ptD_quark[3];
TGraph2D* ptD_gluon[3];
char str[1023];
for(int k=0; k<2;k++)
	{
	nCharged_quark[k]=new TGraph2D();
	nCharged_gluon[k]=new TGraph2D();
	sprintf(str,"nCharged%d_quark",k);
	nCharged_quark[k]->SetName(str);
	sprintf(str,"nCharged%d_gluon",k);
	nCharged_gluon[k]->SetName(str);
	}
for(int k=0; k<2;k++)
	{
	nNeutral_quark[k]=new TGraph2D();
	nNeutral_gluon[k]=new TGraph2D();
	sprintf(str,"nNeutral%d_quark",k);
	nNeutral_quark[k]->SetName(str);
	sprintf(str,"nNeutral%d_gluon",k);
	nNeutral_gluon[k]->SetName(str);
	}
for(int k=0; k<3;k++)
	{
	ptD_quark[k]=new TGraph2D();
	ptD_gluon[k]=new TGraph2D();
	sprintf(str,"ptD%d_quark",k);
	ptD_quark[k]->SetName(str);
	sprintf(str,"ptD%d_gluon",k);
	ptD_gluon[k]->SetName(str);
	}
fprintf(stderr,"Getting Bins\n");
double RhoBins[25];
double PtBins[25];
getBins(PtBins,20,15,1000,true);
getBins(RhoBins,20,0,20,false);
fprintf(stderr,"Starting Loop\n");
for(int j=0; j<3;j++)
{
int count=0;
fprintf(stderr,"VarName=%s\n",VarNames[j].Data());

//printing on the txt file:
if(VarNames[j].Data()[0]!='p') fprintf(fw,"{2 JetPt Rho 1 x [ TMath::Exp( - x *[0]/[1] ) * TMath::Power(x,[0]-1) * TMath::Power([1]/[0],-[0])/TMath::Gamma([0])] QGL %s}\n",VarNames[j].Data()); //TXT
else fprintf(fw,"{2 JetPt Rho 1 x [((x-[0])<0)?0:TMath::Exp( - (x-[0]) *[1]/[2]) * TMath::Power((x-[0]),[1]-1) * TMath::Power([2]/[1],-[1])/TMath::Gamma([1])] QGL %s}\n",VarNames[j].Data()); //TXT 

//il set dei parametri e' messo qui in modo da 'seguire' i risultati dei fit precedenti
gammadistr->SetParLimits(0,1,20);
gammadistr->SetParLimits(1,1,50);

for(int RhoBin=0;RhoBin<20;RhoBin++)
for(int PtBin=0; PtBin<20;PtBin++)
{
	fprintf(stderr,"Bin: Rho %.0lf Pt %.0lf - %.0lf\n",floor(RhoBins[RhoBin]),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]));

//
	functionPtD->SetParameter(0,0.2);
	functionPtD->SetParLimits(0,0.,0.4);//offset

	functionPtD->SetParameter(1,5);//sigma gaus
	functionPtD->SetParLimits(1,1,50);

	functionPtD->SetParameter(2,0.2);
	functionPtD->SetParLimits(2,0.001,0.99);

char plotName[1023];
	sprintf(plotName,"rhoBins_pt%.0f_%.0f/%s_quark_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),VarNames[j].Data(),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin]));
	TH1F* Histo_quark=(TH1F*)f->Get(plotName);	
	sprintf(plotName,"rhoBins_pt%.0f_%.0f/%s_gluon_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),VarNames[j].Data(),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin]));
	TH1F* Histo_gluon=(TH1F*)f->Get(plotName);

//Normalize
Histo_quark->Scale(1./Histo_quark->Integral("width"));
Histo_gluon->Scale(1./Histo_gluon->Integral("width"));

float Max=0;
Max=TMath::Max(Histo_quark->GetMaximum(),Histo_gluon->GetMaximum());
Histo_quark->SetMaximum(Max*1.1);

Histo_quark->SetLineColor(kBlack);
Histo_quark->SetFillColor(kGray);

Histo_gluon->SetLineColor(kRed);
Histo_gluon->SetFillColor(kRed-9);

F->cd();
TCanvas *c1 =new TCanvas();
//c1->SetLogy();
Histo_quark->Draw();
Histo_gluon->Draw("SAME");
f->cd();
TH1F*tmp=(TH1F*)Histo_quark->Clone("tmp"); tmp->SetFillColor(0);tmp->Draw("SAME");
F->cd();
if( (VarNames[j].Data())[0] =='p')
{
	Histo_quark->Fit("functionPtD","N");
	Histo_quark->Fit("functionPtD","NM");	
	functionPtD->SetLineColor(kBlack);
	functionPtD->DrawCopy("SAME");
}
else
	{
	Histo_quark->Fit("gamma","N");//N=Don't Draw
	Histo_quark->Fit("gamma","NM");//N=Don't Draw M=More
	gammadistr->SetLineColor(kBlack);
	gammadistr->DrawCopy("SAME");
	}
//filling TGraph
	switch((VarNames[j].Data())[1])
	{
	case 't':ptD_quark[0]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,functionPtD->GetParameter(0));
		 ptD_quark[1]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,functionPtD->GetParameter(1));
		 ptD_quark[2]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,functionPtD->GetParameter(2));break;
	case 'C':nCharged_quark[0]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(0));
		 nCharged_quark[1]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(1));break;
	case 'N':nNeutral_quark[0]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(0));
		 nNeutral_quark[1]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(1));break;	
	}
	
	if(VarNames[j].Data()[0]!='p')fprintf(fw,"%.0f %.0f %.0f %.0f 4 0 100 %.3f %.3f Q\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],gammadistr->GetParameter(0),gammadistr->GetParameter(1));
	if(VarNames[j].Data()[0]=='p')fprintf(fw,"%.0f %.0f %.0f %.0f 5 0 1 %.3f %.3f %.3f Q\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],functionPtD->GetParameter(0),functionPtD->GetParameter(1),functionPtD->GetParameter(2));
// GLUON 
if( (VarNames[j].Data())[0] =='p')
{
	Histo_gluon->Fit("functionPtD","N");
	Histo_gluon->Fit("functionPtD","NM");
	functionPtD->SetLineColor(kRed);
	functionPtD->DrawCopy("SAME");
}
else
	{
	Histo_gluon->Fit("gamma","N");//N=Don't Draw
	Histo_gluon->Fit("gamma","NM");//N=Don't Draw M=More
	gammadistr->SetLineColor(kRed);
	gammadistr->DrawCopy("SAME");
	}
	switch((VarNames[j].Data())[1])
	{
	case 't':ptD_gluon[0]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,functionPtD->GetParameter(0));
		 ptD_gluon[1]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,functionPtD->GetParameter(1));
		 ptD_gluon[2]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,functionPtD->GetParameter(2));break;
	case 'C':nCharged_gluon[0]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(0));
		 nCharged_gluon[1]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(1));break;
	case 'N':nNeutral_gluon[0]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(0));
		 nNeutral_gluon[1]->SetPoint(count,(PtBins[PtBin]+PtBins[PtBin+1])/2.0,(RhoBins[RhoBin]+RhoBins[RhoBin+1])/2.0,gammadistr->GetParameter(1));break;	
	}
	if(VarNames[j].Data()[0]!='p')fprintf(fw,"%.0f %.0f %.0f %.0f 4 0 100 %.3f %.3f G\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],gammadistr->GetParameter(0),gammadistr->GetParameter(1));
	if(VarNames[j].Data()[0]=='p')fprintf(fw,"%.0f %.0f %.0f %.0f 5 0 1 %.3f %.3f %.3f G\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],functionPtD->GetParameter(0),functionPtD->GetParameter(1),functionPtD->GetParameter(2));




count++;



sprintf(plotName,"%s_pt_%.0lf_%.0lf_rho%.0lf",VarNames[j].Data(),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin]));
c1->SetName(plotName);
c1->Write();
delete c1;

}//for each bin

}//for each var
F->cd();
for(int k=0; k<2;k++)
	{
	nCharged_quark[k]->Write();
	nCharged_gluon[k]->Write();
	nNeutral_quark[k]->Write();
	nNeutral_gluon[k]->Write();
	}
for(int k=0; k<3;k++)
	{
	ptD_quark[k]->Write();
	ptD_gluon[k]->Write();
	}
F->Write();
fprintf(stderr,"Int=%lf\n",functionPtD->Integral(0,1));
return ;
}
