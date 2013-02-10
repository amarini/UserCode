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

int Plot( TH1F* data, vector<TH1F*> &mc, const char * varName,const char*plotDir, int chid2,float lumi){

if(mc.size()<1) return 1;
for(int i=0;i<mc.size();i++) mc[i]->Scale(lumi);
TH1F* varall=(TH1F*)mc[0]->Clone("varall");
for(int i=1;i<mc.size();i++) varall->Add( mc[i] );

if(mc.size()>0) { mc[0]->SetLineColor(kBlack);mc[0]->SetFillColor(kYellow);    mc[0]->SetLineWidth(2);}
if(mc.size()>1) { mc[1]->SetLineColor(kBlack);mc[1]->SetFillColor(kMagenta+2); mc[1]->SetLineWidth(2); }
if(mc.size()>2) { mc[2]->SetLineColor(kBlack);mc[2]->SetFillColor(kGreen);     mc[2]->SetLineWidth(2); }
if(mc.size()>3) { mc[3]->SetLineColor(kBlack);mc[3]->SetFillColor(kGreen+1);   mc[3]->SetLineWidth(2); }
if(mc.size()>4) { mc[4]->SetLineColor(kBlack);mc[4]->SetFillColor(kGreen+2);   mc[4]->SetLineWidth(2); }
if(mc.size()>5) { mc[5]->SetLineColor(kBlack);mc[5]->SetFillColor(kGreen+3);   mc[5]->SetLineWidth(2); }
data->SetMarkerStyle(20);

//mc[0]->GetXaxis()->SetTitle(varName);
//mc[0]->GetYaxis()->SetTitle("events");

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
	
	h->Draw();
	data->Draw("P SAME");
	
	//After it is drawn
	h->SetTitle(mc[0]->GetTitle());
	h->GetXaxis()->SetTitle(mc[0]->GetXaxis()->GetTitle());
	h->GetYaxis()->SetTitle(mc[0]->GetYaxis()->GetTitle());
		
	lowerPad->cd();
	TH1F * R=(TH1F*)data->Clone("ratio");R->Divide(varall);
  	R->GetXaxis()->SetTitle("");
        R->GetYaxis()->SetTitle("Ratio (Data/MC RECO)");
        R->GetYaxis()->SetTitleSize(0.08);
        R->GetYaxis()->SetTitleOffset(0.6);
        R->GetYaxis()->SetLabelSize(0.1);
        R->GetXaxis()->SetLabelSize(0.1);
        R->GetYaxis()->SetRangeUser(0.4,1.6);
        R->SetMarkerStyle(20);
        R->SetMarkerColor(kBlack);
        R->SetMarkerSize(0.7);
        R->SetLineColor(kBlack);
	R->Draw("P");
	lowerPad->SetGridy(2);

c->SaveAs(Form("%s%s_%d.pdf",plotDir,varName,-chid2));
}


int OpenFiles(const char *outDir, const char *plotDir,int chid2=-4,float lumi=1){
        gROOT->SetStyle("Plain");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);

TFile *DY=TFile::Open( Form("%sDY_%d.root",outDir,-chid2) );
TFile *TT=TFile::Open( Form("%sTT_%d.root",outDir,-chid2) );
TFile *WJ=TFile::Open( Form("%sWJ_%d.root",outDir,-chid2) );
TFile *WW=TFile::Open( Form("%sWW_%d.root",outDir,-chid2) );
TFile *WZ=TFile::Open( Form("%sWZ_%d.root",outDir,-chid2) );
TFile *ZZ=TFile::Open( Form("%sZZ_%d.root",outDir,-chid2) );

TFile *DoubleMu=(chid2!=-4)?NULL:TFile::Open( Form("%sDoubleMu_%d.root",outDir,-chid2) );
TFile *DoubleE =(chid2!=-1)?NULL:TFile::Open( Form("%sDoubleE_%d.root",outDir,-chid2) );
TFile *MuEG    =(chid2!=-2)?NULL:TFile::Open( Form("%sMuEG_%d.root",outDir,-chid2) );

vector<string> histoName;
	histoName.push_back("nVtx");
	histoName.push_back("nLeptons");
	histoName.push_back("nPhotons");
	histoName.push_back("nJets");
	histoName.push_back("rho");
	histoName.push_back("rho25");
	histoName.push_back("llM");
	histoName.push_back("llPt");
	histoName.push_back("llY");
	histoName.push_back("llPhi");
	histoName.push_back("lep0Pt");
	histoName.push_back("lep0Eta");
	histoName.push_back("lep0Phi");
	histoName.push_back("lep1Pt");
	histoName.push_back("lep1Eta");
	histoName.push_back("lep1Phi");
	histoName.push_back("jet0Pt");
	histoName.push_back("jet0Eta");
	histoName.push_back("jet0Phi");
	histoName.push_back("jet0QGL");
	histoName.push_back("jet0Btag");
	histoName.push_back("jet1Pt");
	histoName.push_back("jet1Eta");
	histoName.push_back("jet1Phi");
	histoName.push_back("jet1QGL");
	histoName.push_back("jet1Btag");
	histoName.push_back("JetLLDPhi0");
	histoName.push_back("Sum3j");

for(int i=0;i<histoName.size();i++)
	{
	vector<TH1F*> mc;
	TH1F*data;
		switch (chid2)
			{
			case -4:data=(TH1F*)DoubleMu->Get(histoName[i].c_str())->Clone( "data" );break;
			case -2:data=(TH1F*)MuEG->Get(histoName[i].c_str())->Clone( "data" );break;
			case -1:data=(TH1F*)DoubleE->Get(histoName[i].c_str())->Clone( "data" );break;
			}
	mc.push_back( (TH1F*)DY->Get(histoName[i].c_str())->Clone("DY")  );
	mc.push_back( (TH1F*)TT->Get(histoName[i].c_str())->Clone("TT")  );
	mc.push_back( (TH1F*)WJ->Get(histoName[i].c_str())->Clone("WJ")  );
	mc.push_back( (TH1F*)WW->Get(histoName[i].c_str())->Clone("WW")  );
	mc.push_back( (TH1F*)WZ->Get(histoName[i].c_str())->Clone("WZ")  );
	mc.push_back( (TH1F*)ZZ->Get(histoName[i].c_str())->Clone("ZZ")  );

	Plot(data,mc,histoName[i].c_str(),plotDir,chid2,lumi);
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
