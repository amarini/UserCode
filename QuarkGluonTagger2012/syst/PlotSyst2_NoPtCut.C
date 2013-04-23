

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <cstdio>
#include "Syst2.1.C"

#include "TText.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"

using namespace std;
#include "Overflow.h"

int PlotSyst2_NoPtCut(float ptmin=30,float ptmax=50,float rhomin=0, float rhomax=15,float etamin=0, float etamax=2.0,const char*var="QGLHisto")
{

gStyle->SetOptTitle(kFALSE);
gStyle->SetOptStat(kFALSE);

string Zselection="&& axis1_QCJet0>0 && axis2_QCJet0>0 && mZ>70 && mZ<110 && nPFCand_QC_ptCutJet>0 && ptD_QCJet0>0 && (betaStarJet0 < 0.2 * TMath::Log( nvertex - 0.67)) && (deltaPhi_jet>3.1415-0.5) ";

//----- global vars called by minuit :(
TChain *t_mc;
TChain *t_data;

float PtMin=ptmin;
float PtMax=ptmax;
float RhoMin=rhomin;
float RhoMax=rhomax;
float EtaMin=etamin;
float EtaMax=etamax;

string selection=Form("ptJet0>%.0f && ptJet0<%.0f && %.0f<=rhoPF && rhoPF<%.0f && %.1f<=abs(etaJet0) && abs(etaJet0)<%.1f",PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax)+Zselection;
//string varName="QGLikelihood2012Jet0";
//string varName="QGLHisto";

 t_data=new TChain("tree_passedEvents");
 t_mc=new TChain("tree_passedEvents");
	
 t_data->Add("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_Double*.root");
 t_mc->Add("/Users/andreamarini/Documents/QGDiscriminator/ZJet/ZJet_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_2.root");

string sel=string("PUReWeight*eventWeight*("+selection+")");

//INIT
string varName(var);
int nBins=30;
float xMin=0,xMax=1.000001;
if(EtaMax<2.5)
	{t_data->Draw( Form("%s>>h_data(%d,%lf,%lf)",varName.c_str(),nBins,xMin,xMax),selection.c_str(), "E");
	printf("%s Central\n",varName.c_str());}
else
	{t_data->Draw( Form("%sFwd>>h_data(%d,%lf,%lf)",varName.c_str(),nBins,xMin,xMax),selection.c_str(), "E");
	printf("%s Fwd\n",varName.c_str());}
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
float QGLMLP;t_mc->SetBranchAddress("QGLMLP",&QGLMLP);
float eventWeight;t_mc->SetBranchAddress("eventWeight",&eventWeight);
float PUReWeight;t_mc->SetBranchAddress("PUReWeight",&PUReWeight);


TH1F* h_mc1=new TH1F("h_mc1","h_mc1",nBins,xMin,xMax);
TH1F* h_mc0_redone=new TH1F("h_mc0_redone","h_mc0_redone",nBins,xMin,xMax);
h_mc1->Sumw2();
	
TH1F* cut_flow=new TH1F("cut","cut",10,0,10);
for(long  i=0;i<t_mc->GetEntries();i++)
	{
	t_mc->GetEntry(i);	
	//if((i&1023)==1)printf("Entry i=%ld of %ld\n",long(i),long(t_mc->GetEntries()));
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
	//if( !(fabs(ptZ-ptJet0)/(ptZ+ptJet0)<.4 ) ) continue;
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
	
	h_mc1->Fill( Syst(varName.c_str(), ptJet0,  rhoPF,  etaJet0,(varName=="QGLHisto")?QGLHisto:QGLMLP),eventWeight*PUReWeight);
	h_mc0_redone->Fill( (varName=="QGLHisto")?QGLHisto:QGLMLP,eventWeight*PUReWeight);

	}
TCanvas *c0=new TCanvas("c0","c0",800,800);
	bool DrawRatio=true;
	TPad *up=new TPad("up","up",0,0.2,1,1);
	TPad *dn=new TPad("dn","dn",0,0.0,1,.2);
	
	if(DrawRatio){
		up->Draw();
		dn->Draw();
		up->cd();
		}
	
	float data_int=h_data->Integral();
	h_mc0->Scale(data_int*1./h_mc0->Integral());

	h_mc1->Scale(data_int*1./h_mc1->Integral());
	h_data->Scale(data_int*1./h_data->Integral());
	//h_data->Scale(h_mc0->Integral()/h_data->Integral());
	
	h_mc0->GetXaxis()->SetTitle(varName.c_str());
	h_mc0->SetLineColor(kBlue);
	h_mc0->SetLineWidth(2);
	h_mc0->SetFillColor(kBlue-7);
	h_mc0->SetFillStyle(3001);
//	h_mc0_redone->SetLineColor(kGreen);
	h_mc1->SetLineColor(kRed);
	h_mc1->SetLineWidth(2);
	h_data->SetMarkerStyle(20);
	
	double yMin=h_mc0->GetMinimum();
	double yMax=h_mc0->GetMaximum();
	
	yMin=(yMin>h_data->GetMinimum()) ?h_data->GetMinimum() :yMin;
	yMax=(yMax<h_data->GetMaximum()) ?h_data->GetMaximum() :yMax;
	
	yMin=(yMin>h_mc1->GetMinimum()) ?h_mc1->GetMinimum() :yMin;
	yMax=(yMax<h_mc1->GetMaximum()) ?h_mc1->GetMaximum() :yMax;
	
	yMax*=1.1;
	yMin*=0;
	
	h_mc0->SetMaximum(yMax);
	h_mc0->SetMinimum(yMin);

	TH1F* h_mc0_ERRORS=h_mc0->Clone("ERRORS");	
		h_mc0->SetFillStyle(0);
	PaintOverflow(h_mc0_ERRORS,"E2");
	PaintOverflow(h_mc0,"HIST SAME");
//	PaintOverflow(h_mc0_redone,"HIST SAME");
	PaintOverflow(h_mc1,"HIST SAME");
	PaintOverflow(h_data,"P SAME");
	
	TLegend*L=new TLegend(0.4,.5,.6,.8,Form("#splitline{%.0f<P_{T}<%.0f %.0f<#rho<%.0f [GeV]}{%.1f<#eta<%.1f}",PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax) );
	L->AddEntry(h_mc0,"MC","L");
	//L->AddEntry(h_mc0_redone,"MC redone","L");
	L->AddEntry(h_mc1,"MC+syst","L");
	L->AddEntry(h_data,"data","P");
	L->Draw();
		
	if(DrawRatio){
		dn->cd();	
		TH1F* r_mc0=(TH1F*)h_mc0->Clone("r_mc0");r_mc0->Divide(h_data);
		TH1F* r_mc1=(TH1F*)h_mc1->Clone("r_mc1");r_mc1->Divide(h_data);
		r_mc0->GetXaxis()->SetTitleFont(63);r_mc0->GetXaxis()->SetTitleSize(20);
		r_mc0->GetXaxis()->SetTitle("");
		r_mc0->GetYaxis()->SetTitle("Ratio");
		r_mc0->GetYaxis()->SetTitleFont(63);r_mc0->GetYaxis()->SetTitleSize(20);
		r_mc0->GetXaxis()->SetLabelFont(63);r_mc0->GetXaxis()->SetLabelSize(16);
		r_mc0->GetYaxis()->SetLabelFont(63);r_mc0->GetYaxis()->SetLabelSize(16);
		r_mc0->SetMinimum(0.5);
		r_mc0->SetMaximum(1.5);
		PaintOverflow(r_mc0,"P ");
		PaintOverflow(r_mc1,"P SAME");
		TGraph *g=new TGraph();g->SetName("line"); g->SetPoint(0,-1,1);g->SetPoint(1,2,1); g->Draw("L SAME");
		}
	
	c0->SaveAs(Form("Results/%s_noPtCut_%s_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f.pdf",Version.c_str(),varName.c_str(),PtMin,PtMax,RhoMin,RhoMax,EtaMax));


}


