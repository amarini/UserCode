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
bool LogScale=false;
#define MAXBINS 100
float PtBins[MAXBINS];
int nPtBins=0;

const int debug=3;

int Plot( vector<TH1F*> &data, vector<TH1F*> &mc, const char * varName,const char*plotDir, float lumi){
if(debug>1)printf("Start Plot %s\n",varName);

if(mc.size()<1) return 1;
if(data.size()<3) return 1;

for(int i=0;i<mc.size();i++) mc[i]->Scale(lumi);

//Add DY
TH1F* varall=(TH1F*)mc[0]->Clone("varall");
TH1F* data_h=(TH1F*)data[0]->Clone("data");
data_h->Add(data[1]);
data_h->Add(data[2],-1.0);

//for(int i=1;i<mc.size();i++) varall->Add( mc[i] );

if(mc.size()>0) { mc[0]->SetLineColor(kBlack);mc[0]->SetFillColor(kYellow);    mc[0]->SetLineWidth(2);}
if(mc.size()>1) { mc[1]->SetLineColor(kBlack);mc[1]->SetFillColor(kMagenta+2); mc[1]->SetLineWidth(2); }
if(mc.size()>2) { mc[2]->SetLineColor(kBlack);mc[2]->SetFillColor(kGreen);     mc[2]->SetLineWidth(2); }
if(mc.size()>3) { mc[3]->SetLineColor(kBlack);mc[3]->SetFillColor(kGreen+1);   mc[3]->SetLineWidth(2); }
if(mc.size()>4) { mc[4]->SetLineColor(kBlack);mc[4]->SetFillColor(kGreen+2);   mc[4]->SetLineWidth(2); }
if(mc.size()>5) { mc[5]->SetLineColor(kBlack);mc[5]->SetFillColor(kGreen+3);   mc[5]->SetLineWidth(2); }
data_h->SetMarkerStyle(20);

float ymax=TMath::Max(data_h->GetMaximum(),varall->GetMaximum());
mc[0]->SetMaximum(ymax*1.1);

//THStack *h=new THStack("MC_Stack","stack");
//for(int i=mc.size()-1;i>=0;i--) h->Add(mc[i]);
TH1F*h=mc[0];

TCanvas *c=new TCanvas("c1","c1",800,1000);

	TPad*    upperPad = new TPad("upperPad", "upperPad",.005, .25, .995, .995);
        TPad*    lowerPad = new TPad("lowerPad", "lowerPad",.005, .005, .995, .2475);
        upperPad->Draw();
        lowerPad->Draw();
        upperPad->cd();
        if(LogScale)upperPad->SetLogy();

	//printf("going to draw\n");	
	h->Draw("HIST");
	data_h->Draw("P E0 SAME");
	h->Draw("AXIS X+ Y+ SAME");
	h->Draw("AXIS SAME");
	
	//After it is drawn
	h->SetTitle(mc[0]->GetTitle());
	h->GetXaxis()->SetTitle(mc[0]->GetXaxis()->GetTitle());
	h->GetYaxis()->SetTitle(mc[0]->GetYaxis()->GetTitle());

	//printf("going to draw legend\n");	
	TLegend *L=new TLegend(0.15,0.150,0.45,.40,"");
		L->AddEntry(data_h,"data (sf-of)","PF");
		if(mc.size()>0) L->AddEntry(mc[0],"DY","F");
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
	TH1F * R=(TH1F*)data_h->Clone("ratio");
	R->Sumw2();
	R->Divide(varall);
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
	R->Draw("P E0");
	lowerPad->SetGridy(2);

c->SaveAs(Form("%svQG_%s.pdf",plotDir,varName));
}


int OpenFiles(const char *outDir, const char *plotDir,float lumi=1){
if(debug>0)printf("Starting opening files\n");
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(0);

int chid2=-1;
if(debug>0)printf("Opening:%s\n",Form("%svQG_DY_%d.root",outDir,-chid2) );
TFile *DY_e=TFile::Open( Form("%svQG_DY_%d.root",outDir,-chid2) );
TFile *TT_e=TFile::Open( Form("%svQG_TT_%d.root",outDir,-chid2) );
//TFile *WJ_e=TFile::Open( Form("%svQG_WJ_%d.root",outDir,-chid2) );
//TFile *WW_e=TFile::Open( Form("%svQG_WW_%d.root",outDir,-chid2) );
//TFile *WZ_e=TFile::Open( Form("%svQG_WZ_%d.root",outDir,-chid2) );
//TFile *ZZ_e=TFile::Open( Form("%svQG_ZZ_%d.root",outDir,-chid2) );
chid2=-4;
if(debug>0)printf("Opening:%s\n",Form("%svQG_DY_%d.root",outDir,-chid2) );
TFile *DY_mu=TFile::Open( Form("%svQG_DY_%d.root",outDir,-chid2) );
TFile *TT_mu=TFile::Open( Form("%svQG_TT_%d.root",outDir,-chid2) );
//TFile *WJ_mu=TFile::Open( Form("%svQG_WJ_%d.root",outDir,-chid2) );
//TFile *WW_mu=TFile::Open( Form("%svQG_WW_%d.root",outDir,-chid2) );
//TFile *WZ_mu=TFile::Open( Form("%svQG_WZ_%d.root",outDir,-chid2) );
//TFile *ZZ_mu=TFile::Open( Form("%svQG_ZZ_%d.root",outDir,-chid2) );
chid2=-2;
if(debug>0)printf("Opening:%s\n",Form("%svQG_DY_%d.root",outDir,-chid2) );
TFile *DY_emu=TFile::Open( Form("%svQG_DY_%d.root",outDir,-chid2) );
//TFile *TT_emu=TFile::Open( Form("%svQG_TT_%d.root",outDir,-chid2) );
//TFile *WJ_emu=TFile::Open( Form("%svQG_WJ_%d.root",outDir,-chid2) );
//TFile *WW_emu=TFile::Open( Form("%svQG_WW_%d.root",outDir,-chid2) );
//TFile *WZ_emu=TFile::Open( Form("%svQG_WZ_%d.root",outDir,-chid2) );
//TFile *ZZ_emu=TFile::Open( Form("%svQG_ZZ_%d.root",outDir,-chid2) );

TFile *DoubleMu=TFile::Open( Form("%svQG_DoubleMu_%d.root",outDir,4) );
TFile *DoubleE =TFile::Open( Form("%svQG_DoubleE_%d.root",outDir,1) );
TFile *MuEG    =TFile::Open( Form("%svQG_MuEG_%d.root",outDir,2) );

string name;
vector<string> histoName;
if(debug>0)printf("Booking histos\n" );
	//histoName.push_back("nVtx");
	for(int bin=0;bin<nPtBins;bin++)
		{
		name=Form("QGL_jet0_pt%.0f_%.0f",PtBins[bin],PtBins[bin+1])   ;
		histoName.push_back(name);
		name=Form("QGMLP_jet0_pt%.0f_%.0f",PtBins[bin],PtBins[bin+1])   ;
		histoName.push_back(name);
		}

if(debug>0)printf("Histo Loop\n" );
for(int i=0;i<histoName.size();i++)
	{
	vector<TH1F*> mc;
	vector<TH1F*> data;
	
	if(debug>0)printf("Histo data:\n",histoName[i].c_str() );

	data.push_back( (TH1F*)DoubleMu->Get(histoName[i].c_str())->Clone("DoubleMu"));	
	data.push_back( (TH1F*)DoubleE->Get(histoName[i].c_str())->Clone("DoubleE"));	
	data.push_back( (TH1F*)MuEG->Get(histoName[i].c_str())->Clone("MuEG"));	

	if(debug>0)printf("Histo mc DY:\n",histoName[i].c_str() );
	mc.push_back( (TH1F*)DY_mu->Get(histoName[i].c_str())->Clone("DY")  );
	//mc.push_back( (TH1F*)TT_mu->Get(histoName[i].c_str())->Clone("TT")  );
	//mc.push_back( (TH1F*)WJ_mu->Get(histoName[i].c_str())->Clone("WJ")  );
	//mc.push_back( (TH1F*)WW_mu->Get(histoName[i].c_str())->Clone("WW")  );
	//mc.push_back( (TH1F*)WZ_mu->Get(histoName[i].c_str())->Clone("WZ")  );
	//mc.push_back( (TH1F*)ZZ_mu->Get(histoName[i].c_str())->Clone("ZZ")  );

	if(debug>0)printf("Histo mc2:\n",histoName[i].c_str() );
	mc[0]->Add((TH1F*)DY_e->Get(histoName[i].c_str()));
	//mc[1]->Add((TH1F*)TT_e->Get(histoName[i].c_str()));
	//mc[2]->Add((TH1F*)WJ_e->Get(histoName[i].c_str()));
	//mc[3]->Add((TH1F*)WW_e->Get(histoName[i].c_str()));
	//mc[4]->Add((TH1F*)WZ_e->Get(histoName[i].c_str()));
	//mc[5]->Add((TH1F*)ZZ_e->Get(histoName[i].c_str()));
//	if(histoName[i].find("gM")!=string::npos)
//		{
//		data->Rebin(2);for(int i=0;i<mc.size();i++) mc[i]->Rebin(2);
//		data->GetYaxis()->SetTitle("events/20GeV");for(int i=0;i<mc.size();i++) mc[i]->GetYaxis()->SetTitle("events/20GeV");
//		}

	Plot(data,mc,histoName[i].c_str(),plotDir,lumi);
	}

}

#ifdef STANDALONE
#include "src/ReadParameters.C"
int main(int argc,char **argv)
{
//Read A("data/config.ini");
string configFile;
string configFile2;

if(argc<2) configFile="data/config.ini";
else configFile=argv[1];

bool secondConfig=false;
if(argc>=3){configFile2=argv[2];secondConfig=true;}
else configFile2="";

Read A(configFile.c_str());

string DirOut;
string PlotDir;
int chid2; 
float lumi; 

	DirOut=A.ReadParFromMultFile(configFile2.c_str(),"OUTDIR"); DirOut+="/";
	PlotDir=A.ReadParFromMultFile(configFile2.c_str(),"PLOTDIR"); PlotDir+="/";
//	sscanf(A.ReadParFromMultFile(configFile2.c_str(),"CHID2"),"%d",&chid2);
	sscanf(A.ReadParFromMultFile(configFile2.c_str(),"LUMI"),"%f",&lumi);

	const char *ptbins_str=A.ReadParFromMultFile(configFile2.c_str(),"PTBINS");
		nPtBins=0;int n;
		while( sscanf(ptbins_str,"%f%n",&PtBins[nPtBins],&n)>0){ptbins_str+=n;nPtBins++;}
		nPtBins--;

OpenFiles(DirOut.c_str(),PlotDir.c_str(),lumi);
}
#endif 
