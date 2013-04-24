#include "TFile.h"
#include "TH1F.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"
#include "../ZAnalysis/ZAnalysis.h"
#include "TCanvas.h"
#include "math.h"

TH1F * Unfold(TH1F* h1,RooUnfoldResponse *R,int par,int nToys=1000,TH1D **hD=NULL,TMatrixD**Cov=NULL)
{
	RooUnfoldSvd U(R,h1,par,nToys);U.SetNToys(nToys); //toys
	TH1F *hTrueBin = (TH1F*)U.Hreco(RooUnfold::kCovToy);
	if(hD) (*hD)=U.Impl()->GetD();
	if(Cov){
        (*Cov) = (TMatrixD*)U.Ereco(RooUnfold::kCovToy).Clone("cov matrix");
	}
return hTrueBin;
}

int ZUnfolding()
{
//Open Relevant Files
ZAnalysis *A=new ZAnalysis();

string OutputDir=A->Read("data/config.ini","OUTDIR");
float lumi;
sscanf(A->Read("data/config.ini","LUMI"),"%f",&lumi);

//OPEN OUTPUT
TFile *Out=TFile::Open(Form("%s/unfold.root",OutputDir.c_str()),"RECREATE");

// DATA
TFile *DY_4=TFile::Open( (OutputDir+"/output_type1_chid4.root").c_str());
TFile *DY_1=TFile::Open( (OutputDir+"/output_type1_chid1.root").c_str());

TFile *TT_4=TFile::Open( (OutputDir+"/output_type2_chid4.root").c_str());
TFile *TT_1=TFile::Open( (OutputDir+"/output_type2_chid1.root").c_str());

TFile *WW_4=TFile::Open( (OutputDir+"/output_type6_chid4.root").c_str());
TFile *WW_1=TFile::Open( (OutputDir+"/output_type6_chid1.root").c_str());

TFile *WZ_4=TFile::Open( (OutputDir+"/output_type4_chid4.root").c_str());
TFile *WZ_1=TFile::Open( (OutputDir+"/output_type4_chid1.root").c_str());

TFile *WJ_4=TFile::Open( (OutputDir+"/output_type3_chid4.root").c_str());
TFile *WJ_1=TFile::Open( (OutputDir+"/output_type3_chid1.root").c_str());

TFile *ZZ_4=TFile::Open( (OutputDir+"/output_type5_chid4.root").c_str());
TFile *ZZ_1=TFile::Open( (OutputDir+"/output_type5_chid1.root").c_str());

TFile *DoubleMu=TFile::Open( (OutputDir+"/output_type0_chid4.root").c_str());
TFile *DoubleE =TFile::Open( (OutputDir+"/output_type0_chid1.root").c_str());
TFile *EMu     =TFile::Open( (OutputDir+"/output_type0_chid2.root").c_str());

//Get Response Matrixes
//Use this code if you have the response directly: they cannot be merged
//RooUnfoldResponse *R_4=(RooUnfoldResponse*)DY_4->Get("R_llPt")->Clone("R_llPt_4");
//RooUnfoldResponse *R_1=(RooUnfoldResponse*)DY_4->Get("R_llPt")->Clone("R_llPt_1");
TH1F* R_4_GEN= (TH1F*)DY_4->Get("R_llPt_GEN")->Clone("R_llPt_GEN_4");
TH1F* R_4_RECO= (TH1F*)DY_4->Get("R_llPt_RECO")->Clone("R_llPt_RECO_4");
TH2F* R_4_MATRIX= (TH2F*)DY_4->Get("R_llPt_MATRIX")->Clone("R_llPt_MATRIX_4");

TH1F* R_1_GEN= (TH1F*)DY_1->Get("R_llPt_GEN")->Clone("R_llPt_GEN_1");
TH1F* R_1_RECO= (TH1F*)DY_1->Get("R_llPt_RECO")->Clone("R_llPt_RECO_1");
TH2F* R_1_MATRIX= (TH2F*)DY_1->Get("R_llPt_MATRIX")->Clone("R_llPt_MATRIX_1");

RooUnfoldResponse *R_4=new RooUnfoldResponse(R_4_RECO,R_4_GEN,R_4_MATRIX);
RooUnfoldResponse *R_1=new RooUnfoldResponse(R_1_RECO,R_1_GEN,R_1_MATRIX);
//Get Histos

TH1F*h_mu=(TH1F*)DoubleMu->Get("llPt")->Clone("llPt_4");
TH1F*h_e =(TH1F*)DoubleE ->Get("llPt")->Clone("llPt_1");
TH1F*h_emu=(TH1F*)EMu    ->Get("llPt")->Clone("llPt_2");

//JES UP/DN
TH1F*h_mu_JES_UP =(TH1F*)DoubleMu->Get("llPt_JES_UP")->Clone("llPt_JES_UP_4");
TH1F*h_e_JES_UP  =(TH1F*)DoubleE ->Get("llPt_JES_UP")->Clone("llPt_JES_UP_1");
TH1F*h_emu_JES_UP=(TH1F*)EMu     ->Get("llPt_JES_UP")->Clone("llPt_JES_UP_2");
TH1F*h_mu_JES_DN =(TH1F*)DoubleMu->Get("llPt_JES_DN")->Clone("llPt_JES_DN_4");
TH1F*h_e_JES_DN  =(TH1F*)DoubleE ->Get("llPt_JES_DN")->Clone("llPt_JES_DN_1");
TH1F*h_emu_JES_DN=(TH1F*)EMu     ->Get("llPt_JES_DN")->Clone("llPt_JES_DN_2");

//MC
TH1F*h_dy_4=(TH1F*)DY_4->Get("llPt")->Clone("llPt_dy_4");
TH1F*h_dy_1=(TH1F*)DY_1->Get("llPt")->Clone("llPt_dy_1");

TH1F*h_tt_4=(TH1F*)TT_4->Get("llPt")->Clone("llPt_tt_4");
TH1F*h_tt_1=(TH1F*)TT_1->Get("llPt")->Clone("llPt_tt_1");

TH1F*h_vv_4=(TH1F*)WJ_4->Get("llPt")->Clone("llPt_vv_4");
TH1F*h_vv_1=(TH1F*)WJ_1->Get("llPt")->Clone("llPt_vv_1");
	h_vv_4->Add((TH1F*)WW_4->Get("llPt"));
	h_vv_1->Add((TH1F*)WW_1->Get("llPt"));

	h_vv_4->Add((TH1F*)WZ_4->Get("llPt"));
	h_vv_1->Add((TH1F*)WZ_1->Get("llPt"));

	h_vv_4->Add((TH1F*)ZZ_4->Get("llPt"));
	h_vv_1->Add((TH1F*)ZZ_1->Get("llPt"));

TH1F* h_dy_gen=(TH1F*)DY_4->Get("llPtGEN")->Clone("llPt_dy_gen");
	h_dy_gen->Add((TH1F*)DY_1->Get("llPtGEN"));

//------------------LUMI--------------------------------
h_dy_4->Scale(lumi);
h_dy_1->Scale(lumi);
h_tt_4->Scale(lumi);
h_tt_1->Scale(lumi);
h_vv_4->Scale(lumi);
h_vv_1->Scale(lumi);
h_dy_gen->Scale(lumi);
//------------------------------------------------------

//----------------BKG---SUB-----------------------------
h_mu->Add(h_emu,-0.5);
h_e ->Add(h_emu,-0.5);
h_mu->Add(h_vv_4,-1.0);
h_e ->Add(h_vv_1,-1.0);

h_mu_JES_UP->Add(h_emu_JES_UP,-0.5);
h_e_JES_UP ->Add(h_emu_JES_UP,-0.5);
h_mu_JES_UP->Add(h_vv_4,-1.0);
h_e_JES_UP ->Add(h_vv_1,-1.0);

h_mu_JES_DN->Add(h_emu_JES_DN,-0.5);
h_e_JES_DN ->Add(h_emu_JES_DN,-0.5);
h_mu_JES_DN->Add(h_vv_4,-1.0);
h_e_JES_DN ->Add(h_vv_1,-1.0);
//----------------UNFOLD--------------------------------------
int ParUnf=30;
TH1F *u_mu=Unfold(h_mu,R_4,ParUnf);u_mu->SetName("u_mu");
TH1F *u_e=Unfold(h_e,R_1,ParUnf);u_e->SetName("u_e");
//----------------UNFOLD STUDIES PAR --------------------
	map<string,TH1F*> unfold_pars;
	Out->cd();
	TCanvas *c=new TCanvas("Par_studies","par_studies",800,600);
	for(int p=10;p<h_mu->GetNbinsX()-1;p+=3){
	  string name=Form("mu_%d",p);
	  unfold_pars[name]=Unfold(h_mu,R_4,p,100);
	  unfold_pars[name]->SetName(name.c_str());
	  unfold_pars[name]->SetLineColor(kGreen+2*p);	
	  if(p!=1) unfold_pars[name]->Draw("P E SAME");	
		else unfold_pars[name]->Draw("P E ");
	
	  name=Form("e_%d",p);
	  unfold_pars[name]=Unfold(h_e,R_1,p,100);
	  unfold_pars[name]->SetName(name.c_str());
	  unfold_pars[name]->SetLineColor(kGreen+2*p+1);	
	  unfold_pars[name]->Draw("P E SAME");	
	}
	Out->cd();
	c->Write();
//---------------- END UNFOLD STUDIES --------------------

//-------------------SYST-----------------------------------
TH1F *u_mu_JES_UP=Unfold(h_mu_JES_UP,R_4,ParUnf); u_mu_JES_UP->SetName("u_mu_JES_UP");
TH1F *u_e_JES_UP =Unfold(h_e_JES_UP,R_4,ParUnf);  u_e_JES_UP->SetName("u_e_JES_UP");
TH1F *u_mu_JES_DN=Unfold(h_mu_JES_DN,R_4,ParUnf);u_mu_JES_DN->SetName("u_mu_JES_DN");
TH1F *u_e_JES_DN =Unfold(h_e_JES_DN,R_4,ParUnf); u_e_JES_DN->SetName("u_e_JES_DN");

//-------------------COMBINATION----------------------------
TH1F*u_l=(TH1F*)u_mu->Clone("llPt_l");
	u_l->Add(u_e);
TH1F*u_l_JES_UP=(TH1F*)u_mu_JES_UP->Clone("llPt_JES_UP_l");
	u_l_JES_UP->Add(u_e_JES_UP);
TH1F*u_l_JES_DN=(TH1F*)u_mu_JES_DN->Clone("llPt_JES_DN_l");
	u_l_JES_DN->Add(u_e_JES_DN);
TH1F*syst=(TH1F*)u_l->Clone("syst");
	for(int i=0;i<=syst->GetNbinsX()+1;i++)
		{
		syst->SetBinContent(i,     (u_l_JES_UP->GetBinContent(i)+u_l_JES_DN->GetBinContent(i) )/2.0);
		syst->SetBinError  (i, fabs(u_l_JES_UP->GetBinContent(i)-u_l_JES_DN->GetBinContent(i) )/2.0);
		}
//-------------------CANVAS---------------------------------
TCanvas *c1=new TCanvas("canvas","canvas",800,800);
h_dy_gen->SetLineColor(kBlue);h_dy_gen->SetLineWidth(2);
u_l->SetMarkerStyle(20);
syst->SetFillColor(kRed);syst->SetFillStyle(3001); syst->SetMarkerStyle(0);

u_l->Draw("AXIS");
u_l->Draw("AXIS X+ Y+ SAME");
syst->Draw("E4 SAME");
h_dy_gen->Draw("HIST SAME");
u_l->Draw("P SAME");

Out->cd();
c1->Write();
h_mu->Write();
h_e->Write();
h_dy_gen->Write();
u_mu->Write();
u_mu_JES_UP->Write();
u_mu_JES_DN->Write();
u_e->Write();
u_e_JES_UP->Write();
u_e_JES_DN->Write();
u_l->Write();
u_l_JES_UP->Write();
u_l_JES_DN->Write();
syst->Write();

//write everything
//..
}


int main(){
ZUnfolding();
}
