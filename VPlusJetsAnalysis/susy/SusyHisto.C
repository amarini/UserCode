#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

using namespace std;
bool LogScale=true;

int Plot( TH1F* data, TH1F*data2, vector<TH1F*> &mc, const char * varName,const char*plotDir,float lumi,const char *LegTitle=""){

if(mc.size()<1) return 1;
for(int i=0;i<mc.size();i++) mc[i]->Scale(lumi);
TH1F* varall=(TH1F*)mc[0]->Clone("varall");
for(int i=1;i<mc.size();i++) varall->Add( mc[i] );



int bin0=varall->FindBin(91-5);
int bin1=varall->FindBin(91+5);
for(int i=0;i<mc.size();i++) mc[i]->Scale(data->Integral(bin0,bin1)/varall->Integral(bin0,bin1));
varall->Scale(data->Integral(bin0,bin1)/varall->Integral(bin0,bin1));
mc[0]->GetYaxis()->SetTitle("events Norm. to the Z peak");

if(mc.size()>0) { mc[0]->SetLineColor(kBlack);mc[0]->SetFillColor(kYellow);    mc[0]->SetLineWidth(2);}
if(mc.size()>1) { mc[1]->SetLineColor(kBlack);mc[1]->SetFillColor(kMagenta+2); mc[1]->SetLineWidth(2); }
if(mc.size()>2) { mc[2]->SetLineColor(kBlack);mc[2]->SetFillColor(kGreen);     mc[2]->SetLineWidth(2); }
if(mc.size()>3) { mc[3]->SetLineColor(kBlack);mc[3]->SetFillColor(kGreen+1);   mc[3]->SetLineWidth(2); }
if(mc.size()>4) { mc[4]->SetLineColor(kBlack);mc[4]->SetFillColor(kGreen+2);   mc[4]->SetLineWidth(2); }
if(mc.size()>5) { mc[5]->SetLineColor(kBlack);mc[5]->SetFillColor(kGreen+3);   mc[5]->SetLineWidth(2); }
data->SetMarkerStyle(20);
data2->SetMarkerStyle(24); data2->SetMarkerColor(kRed+2);  data2->SetLineColor(kRed+2);

//mc[0]->GetXaxis()->SetTitle(varName);

float ymax=TMath::Max(data->GetMaximum(),varall->GetMaximum());
mc[0]->SetMaximum(ymax*1.1);

THStack *h=new THStack("MC_Stack","stack");
for(int i=mc.size()-1;i>=0;i--) h->Add(mc[i]);


TCanvas *c=new TCanvas("c1","c1",800,1000);

	TPad*    upperPad = new TPad("upperPad", "upperPad",.005, .25, .995, .995);
        TPad*    lowerPad = new TPad("lowerPad", "lowerPad",.005, .005, .995, .2475);
        upperPad->Draw();
        lowerPad->Draw();
        upperPad->cd();
        if(LogScale)upperPad->SetLogy();

	//printf("going to draw\n");	
	h->Draw("HIST");
	data->Draw("P E0 SAME");
	data2->Draw("P E0 SAME");
	h->Draw("AXIS X+ Y+ SAME");
	h->Draw("AXIS SAME");
	
	//After it is drawn
	h->SetTitle(mc[0]->GetTitle());
	h->GetXaxis()->SetTitle(mc[0]->GetXaxis()->GetTitle());
	h->GetYaxis()->SetTitle(mc[0]->GetYaxis()->GetTitle());

	//printf("going to draw legend\n");	
	TLegend *L=new TLegend(0.75,0.45,0.89,.89,LegTitle);
		L->AddEntry(data,"data SF","PF");
		L->AddEntry(data2,"data OF","PF");
		if(mc.size()>0) L->AddEntry(mc[0],"DY","F");
		if(mc.size()>1) L->AddEntry(mc[1],"TT","F");
		if(mc.size()>2) L->AddEntry(mc[2],"WJ","F");
		if(mc.size()>3) L->AddEntry(mc[3],"WW","F");
		if(mc.size()>4) L->AddEntry(mc[4],"WZ","F");
		if(mc.size()>5) L->AddEntry(mc[5],"ZZ","F");
		L->SetBorderSize(0);
		L->SetFillColor(0);
		L->SetFillStyle(0);
	L->Draw();
	//printf("going to draw latex\n");	
	TLatex *latex=new TLatex();
		latex->SetNDC();
		latex->SetTextFont(63);
		latex->SetTextSize(28);
		latex->SetTextAlign(21);
		latex->DrawLatex(0.5,.91,"CMS, Work in Progress, #sqrt{s}=8 TeV,L=18.7fb^{-1}");
		
	lowerPad->cd();
	TH1F * R=(TH1F*)data->Clone("ratio");
	R->Sumw2();
	R->Divide(varall);
  	R->GetXaxis()->SetTitle("");
        R->GetYaxis()->SetTitle("Ratio (Data/MC RECO)");
        R->GetYaxis()->SetTitleSize(0.08);
        R->GetYaxis()->SetTitleOffset(0.6);
        R->GetYaxis()->SetLabelSize(0.1);
        R->GetXaxis()->SetLabelSize(0.1);
        R->GetYaxis()->SetRangeUser(0.1,5);
        R->SetMarkerStyle(20);
        R->SetMarkerColor(kBlack);
        R->SetMarkerSize(0.7);
        R->SetLineColor(kBlack);
	R->Draw("P E0");
	lowerPad->SetGridy(2);
	TH1F * R2=(TH1F*)data->Clone("ratio");R2->Sumw2(); 
		R2->Divide(data2);
		R2->SetMarkerStyle(24);
        	R2->SetMarkerSize(0.7);
        	R2->SetMarkerColor(kRed+2);
        	R2->SetLineColor(kRed+2);
		R2->Draw("P E0 SAME");
	

c->SaveAs(Form("%sSusy_%s.pdf",plotDir,varName));
}


int OpenFiles(const char *outDir, const char *plotDir,int chid2=-4,float lumi=1){
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);

TFile *DY_mu=TFile::Open( Form("%sSusy_DY_%d.root",outDir,4) );
TFile *TT_mu=TFile::Open( Form("%sSusy_TT_%d.root",outDir,4) );
TFile *WJ_mu=TFile::Open( Form("%sSusy_WJ_%d.root",outDir,4) );
TFile *WW_mu=TFile::Open( Form("%sSusy_WW_%d.root",outDir,4) );
TFile *WZ_mu=TFile::Open( Form("%sSusy_WZ_%d.root",outDir,4) );
TFile *ZZ_mu=TFile::Open( Form("%sSusy_ZZ_%d.root",outDir,4) );

TFile *DY_e=TFile::Open( Form("%sSusy_DY_%d.root",outDir,1) );
TFile *TT_e=TFile::Open( Form("%sSusy_TT_%d.root",outDir,1) );
TFile *WJ_e=TFile::Open( Form("%sSusy_WJ_%d.root",outDir,1) );
TFile *WW_e=TFile::Open( Form("%sSusy_WW_%d.root",outDir,1) );
TFile *WZ_e=TFile::Open( Form("%sSusy_WZ_%d.root",outDir,1) );
TFile *ZZ_e=TFile::Open( Form("%sSusy_ZZ_%d.root",outDir,1) );

TFile *DoubleMu=TFile::Open( Form("%sSusy_DoubleMu_%d.root",outDir,4) );
TFile *DoubleE =TFile::Open( Form("%sSusy_DoubleE_%d.root",outDir,1) );
TFile *MuEG    =TFile::Open( Form("%sSusy_MuEG_%d.root",outDir,2) );

vector<string> histoName;
	histoName.push_back("llM");
	histoName.push_back("llM_BS");
	histoName.push_back("llM_2j");
	histoName.push_back("llM_2j_BS");

for(int i=0;i<histoName.size();i++)
	{
	vector<TH1F*> mc;
	TH1F*data0;
	TH1F*data1;
	TH1F*data2;
			data0=(TH1F*)DoubleMu->Get(histoName[i].c_str())->Clone( "data0" );
			data1=(TH1F*)DoubleE->Get(histoName[i].c_str())->Clone( "data1" );
			data2=(TH1F*)MuEG->Get(histoName[i].c_str())->Clone( "data2" );

TH1F* h_dy=(TH1F*)DY_mu->Get(histoName[i].c_str())->Clone("DY");h_dy->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
TH1F* h_tt=(TH1F*)TT_mu->Get(histoName[i].c_str())->Clone("TT");h_tt->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
TH1F* h_wj=(TH1F*)WJ_mu->Get(histoName[i].c_str())->Clone("WJ");h_wj->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
TH1F* h_ww=(TH1F*)WW_mu->Get(histoName[i].c_str())->Clone("WW");h_ww->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
TH1F* h_wz=(TH1F*)WZ_mu->Get(histoName[i].c_str())->Clone("WZ");h_wz->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
TH1F* h_zz=(TH1F*)ZZ_mu->Get(histoName[i].c_str())->Clone("ZZ");h_zz->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));

//TH1F* h_dy=(TH1F*)DY_mu->Get(histoName[i].c_str())->Clone("DY");//h_dy->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
//TH1F* h_tt=(TH1F*)TT_mu->Get(histoName[i].c_str())->Clone("TT");//h_tt->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
//TH1F* h_wj=(TH1F*)WJ_mu->Get(histoName[i].c_str())->Clone("WJ");//h_wj->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
//TH1F* h_ww=(TH1F*)WW_mu->Get(histoName[i].c_str())->Clone("WW");//h_ww->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
//TH1F* h_wz=(TH1F*)WZ_mu->Get(histoName[i].c_str())->Clone("WZ");//h_wz->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
//TH1F* h_zz=(TH1F*)ZZ_mu->Get(histoName[i].c_str())->Clone("ZZ");//h_zz->Add( (TH1F*)DY_e->Get(histoName[i].c_str()));
	
	data0->Rebin(4);
	data1->Rebin(4);
	data2->Rebin(4);
	h_dy ->Rebin(4);
	h_tt ->Rebin(4);
	h_wj ->Rebin(4);
	h_ww ->Rebin(4);
	h_wz ->Rebin(4);
	h_zz ->Rebin(4);
	
	mc.push_back( h_dy  );
	mc.push_back( h_tt  );
	mc.push_back( h_wj  );
	mc.push_back( h_ww  );
	mc.push_back( h_wz  );
	mc.push_back( h_zz  );

	TH1F* data=(TH1F*)data0->Clone("data");	
	data->Add(data1);
	//data2->Scale(2); //SF: e+e- mu+mu- OF e+mu- mu+e-
	Plot(data,data2,mc,histoName[i].c_str(),plotDir,lumi,"e^{+}e^{-}+#mu^{+}#mu{-} vs e^{#pm}#mu^{#pm}");

	data2->Scale(.5);
	Plot(data0,data2,mc,(histoName[i]+"_mu").c_str(),plotDir,lumi,"#mu^{+}#mu{-} vs .5 e^{#pm}#mu^{#pm}");
	Plot(data1,data2,mc,(histoName[i]+"_e").c_str(),plotDir,lumi,"e^{+}e^{-} vs .5 e^{#pm}#mu^{#pm}");
	}

}

#ifdef STANDALONE
#include "src/ReadParameters.C"
int main(int argc,char **argv)
{
//Read A("data/config.ini");
string configFile;

if(argc<2) configFile="data/config.ini";
else configFile=argv[1];

Read A(configFile.c_str());

string DirOut=A.ReadParameter("OUTDIR"); DirOut+="/";
string PlotDir=A.ReadParameter("PLOTDIR"); PlotDir+="/";
int chid2; sscanf(A.ReadParameter("CHID2"),"%d",&chid2);
float lumi; sscanf(A.ReadParameter("LUMI"),"%f",&lumi);

OpenFiles(DirOut.c_str(),PlotDir.c_str(),chid2,lumi);
}
#endif 
