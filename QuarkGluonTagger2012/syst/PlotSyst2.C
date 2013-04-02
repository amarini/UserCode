

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <cstdio>
#include "Syst2.C"

#include "TText.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

using namespace std;
void PaintOverflow(TH1F*h,const char *opt="")
{
   // This function paint the histogram h with an extra bin for overflows

   const char* name  = h->GetName();
   const char* title = h->GetTitle();
   Int_t nx    = h->GetNbinsX()+2;
   Double_t bw = h->GetBinWidth(nx);
   Double_t x1 = h->GetBinLowEdge(1)-bw;
   Double_t x2 = h->GetBinLowEdge(nx)+bw;

   // Book a temporary histogram having ab extra bin for overflows
   TH1F *htmp = new TH1F(Form("%s_overflow",name), title, nx, x1, x2);
	htmp->SetLineColor(h->GetLineColor());
	htmp->SetLineStyle(h->GetLineStyle());
	htmp->SetLineWidth(h->GetLineWidth());
	htmp->SetMarkerColor(h->GetMarkerColor());
	htmp->SetMarkerStyle(h->GetMarkerStyle());

	htmp->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
   // Fill the new hitogram including the extra bin for overflows
   for (Int_t i=1; i<=nx; i++) {
      //htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
	htmp->SetBinContent(i,h->GetBinContent(i-1));
	htmp->SetBinError  (i,h->GetBinError(i-1));
   }

   // Restore the number of entries
   htmp->SetEntries(h->GetEntries());

   // Draw the temporary histogram
   htmp->Draw(opt);
   if(string(opt).find("SAME") ==string::npos)
   {
   TText *t = new TText(x2-bw/2,h->GetBinContent(nx),"Overflow");
   t->SetTextAngle(90);
   t->SetTextAlign(12);
   t->SetTextSize(0.03);;
   t->Draw(); 
	t = new TText(x1+bw/2,h->GetBinContent(nx),"Underflow");
   t->SetTextAngle(90);
   t->SetTextAlign(12);
   t->SetTextSize(0.03);;
   t->Draw(); 
   }
}

int PlotSyst2()
{

gStyle->SetOptTitle(kFALSE);
gStyle->SetOptStat(kFALSE);

string Zselection="&& axis1_QCJet0>0 && axis2_QCJet0>0 && mZ>70 && mZ<110 && abs(ptZ-ptJet0)/(ptZ+ptJet0)<.4 && nPFCand_QC_ptCutJet>0 && ptD_QCJet0>0 && (betaStarJet0 < 0.2 * TMath::Log( nvertex - 0.67)) && (deltaPhi_jet>3.1415-0.5) ";

//----- global vars called by minuit :(
TTree *t_mc;
TTree *t_data;

float PtMin=80;
float PtMax=120;
float RhoMin=0;
float RhoMax=40.;
float EtaMin=3.0;
float EtaMax=4.7;

string selection=Form("ptJet0>%.0f && ptJet0<%.0f && %.0f<=rhoPF && rhoPF<%.0f && %.1f<=abs(etaJet0) && abs(etaJet0)<%.1f",PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax)+Zselection;
//string varName="QGLikelihood2012Jet0";
//string varName="QGLHisto";

 TFile *f_data=TFile::Open("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DoubleMu-Run2012C.root");
 t_data=(TTree*)f_data->Get("tree_passedEvents");
 TFile *f_mc=TFile::Open("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_2.root");
 t_mc=(TTree*)f_mc->Get("tree_passedEvents");

string sel=string("PUReWeight*eventWeight*("+selection+")");

//INIT
string varName("QGLHisto");
int nBins=30;
float xMin=0,xMax=1.000001;
t_data->Draw( Form("%s>>h_data(%d,%lf,%lf)",varName.c_str(),nBins,xMin,xMax),selection.c_str(), "E");
t_mc->Draw( Form("%s>>h_mc0(%d,%lf,%lf)",varName.c_str(),nBins,xMin,xMax),sel.c_str(), "E");

TH1F* h_mc0=(TH1F*)gDirectory->Get("h_mc0");
TH1F* h_data=(TH1F*)gDirectory->Get("h_data");


float axis1_QCJet0;t_mc->SetBranchAddress("axis1_QCJet0",&axis1_QCJet0);
float axis2_QCJet0;t_mc->SetBranchAddress("axis2_QCJet0",&axis2_QCJet0);
int   nPFCand_QC_ptCutJet;t_mc->SetBranchAddress("nPFCand_QC_ptCutJet",&nPFCand_QC_ptCutJet);
float ptD_QCJet0;t_mc->SetBranchAddress("ptD_QCJet0",&ptD_QCJet0);
float betaStarJet0;t_mc->SetBranchAddress("betaStarJet0",&betaStarJet0);
int   nvertex;t_mc->SetBranchAddress("nvertex",&nvertex);
float deltaPhi_jet;t_mc->SetBranchAddress("deltaPhi_jet",&deltaPhi_jet);
float mZ;t_mc->SetBranchAddress("mZ",&mZ);
float ptZ;t_mc->SetBranchAddress("ptZ",&ptZ);
float ptJet0;t_mc->SetBranchAddress("ptJet0",&ptJet0);
float rhoPF;t_mc->SetBranchAddress("rhoPF",&rhoPF);
float etaJet0;t_mc->SetBranchAddress("etaJet0",&etaJet0);
float QGLHisto;t_mc->SetBranchAddress("QGLHisto",&QGLHisto);
float QGLHisto;t_mc->SetBranchAddress("QGLHisto",&QGLHisto);
float eventWeight;t_mc->SetBranchAddress("eventWeight",&eventWeight);
float PUReWeight;t_mc->SetBranchAddress("PUReWeight",&PUReWeight);


TH1F* h_mc1=new TH1F("h_mc1","h_mc1",nBins,xMin,xMax);
TH1F* h_mc0_redone=new TH1F("h_mc0_redone","h_mc0_redone",nBins,xMin,xMax);
h_mc1->Sumw2();
	
TH1F* cut_flow=new TH1F("cut","cut",10,0,10);
for(long long i=0;i<t_mc->GetEntries();i++)
	{
	t_mc->GetEntry(i);	
//string Zselection="&& axis1_QCJet0>0 && axis2_QCJet0>0 && mZ>70 && mZ<110 && abs(ptZ-ptJet0)/(ptZ+ptJet0)<.4 && nPFCand_QC_ptCutJet>0 && ptD_QCJet0>0 && (betaStarJet0 < 0.2 * TMath::Log( nvertex - 0.67)) && (deltaPhi_jet>3.1415-0.5) ";
//	cut_flow->Fill(0);
	if( !(axis1_QCJet0>0) ) continue;
//		cout<<"1 " ;cut_flow->Fill(1);
	if( !(axis2_QCJet0>0) ) continue;
//		cout<<"2 " ;cut_flow->Fill(2);
	if( !(mZ>70) ) continue;
//		cout<<"3 " ;cut_flow->Fill(3);
	if( !(mZ<110) ) continue;
//		cout<<"4 " ;cut_flow->Fill(4);
	if( !(fabs(ptZ-ptJet0)/(ptZ+ptJet0)<.4 ) ) continue;
//		cout<<"5 " ;cut_flow->Fill(5);
	if( !( nPFCand_QC_ptCutJet>0 && ptD_QCJet0>0) ) continue;
//		cout<<"6 " ;cut_flow->Fill(6);
	if( !(betaStarJet0 < 0.2 * TMath::Log( nvertex - 0.67)) ) continue;
//		cout<<"7 " ;cut_flow->Fill(7);
	if( !(deltaPhi_jet>3.1415-0.5)  ) continue;
//		cout<<"8 ";cut_flow->Fill(8);
	//normal selection
	if( !(ptJet0>PtMin && ptJet0<PtMax && RhoMin<=rhoPF && rhoPF<RhoMax && EtaMin<=fabs(etaJet0) && fabs(etaJet0)<EtaMax) ) continue;
//		cout<<"--- "<<endl;cut_flow->Fill(9);
	
	h_mc1->Fill( Syst("QGLHisto", ptJet0,  rhoPF,  etaJet0,QGLHisto),eventWeight*PUReWeight);
	h_mc0_redone->Fill( QGLHisto,eventWeight*PUReWeight);

	}
TCanvas *c0=new TCanvas("c0","c0",800,600);
	h_mc0->Scale(1./h_mc0->Integral());
	h_mc1->Scale(1./h_mc1->Integral());
	h_data->Scale(1./h_data->Integral());
	//h_data->Scale(h_mc0->Integral()/h_data->Integral());
	
	h_mc0->GetXaxis()->SetTitle(varName.c_str());
	h_mc0->SetLineColor(kRed);
//	h_mc0_redone->SetLineColor(kGreen);
	h_mc1->SetLineColor(kBlue);
	h_data->SetMarkerStyle(20);
	
//	h_mc0->Draw("HIST");
//	h_mc1->Draw("HIST SAME");
//	h_data->Draw("P SAME");
	
	PaintOverflow(h_mc0,"HIST");
//	PaintOverflow(h_mc0_redone,"HIST SAME");
	PaintOverflow(h_mc1,"HIST SAME");
	PaintOverflow(h_data,"P SAME");
	
	TLegend*L=new TLegend(0.4,.5,.6,.8,Form("#splitline{%.0f<P_{T}<%.0f %.0f<#rho<%.0f [GeV]}{%.1f<#eta<%.1f}",PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax) );
	L->AddEntry(h_mc0,"MC","L");
	//L->AddEntry(h_mc0_redone,"MC redone","L");
	L->AddEntry(h_mc1,"MC+syst","L");
	L->AddEntry(h_data,"data","P");
	L->Draw();
	
	c0->SaveAs(Form("Results/Syst2_QGLHisto_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f.pdf",PtMin,PtMax,RhoMin,RhoMax,EtaMax));

//TCanvas *c1=new TCanvas("c1","c1",800,600);
//	TH1F* h_mc_mean=(TH1F*)h_mc0->Clone("h_mc_mean");
//	h_mc_mean->Add(h_mc1);
//	h_mc_mean->Scale(0.5);
//	for(int j=0;j<=h_mc_mean->GetNbinsX()+1;j++) h_mc_mean->SetBinError(j,fabs(h_mc0->GetBinContent(j)-h_mc1->GetBinContent(j) )/2.0) ;
//	h_mc_mean->SetFillStyle(3003);
//	h_mc_mean->SetFillColor(kRed);
//	h_mc_mean->SetMarkerStyle(0);
//	h_mc_mean->SetMarkerColor(0);
//	h_mc_mean->Draw("P E4");
//	h_data->Draw("P SAME");

//new TCanvas();
//cut_flow->Draw();

}


