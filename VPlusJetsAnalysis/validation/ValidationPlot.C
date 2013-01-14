#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <string>
#include <map>
#include <vector>
#include <iostream>

using namespace std;
const int debug=1;//0 NO DEBUG 1 MINIMAL 2 MAX

void ValidationPlotData(const char *fileIn,const char *fileOut, const char *dirOut);
int CreateHisto(map<string,TH1F*> &histos);
void ValidationPlot();
void Loop(map<string,TH1F*> &histos, TTree *t,int type=0); //0 data 1 mc



// ------------------ Varibales plots -------------------------
int CreateHisto(map<string,TH1F*> &histos)
{
//DEBUG LEVEL 2
//--- counting
histos["nVtx"] 		= new TH1F("nVtx","nVtx;N_{Vtx};events",30,0,30);
histos["nLeptons"] 	= new TH1F("nLeptons","nLeptons;N_{lep};events",10,0,10);
histos["nPhotons"] 	= new TH1F("nPhotons","nPhotons;N_{#gamma};events",10,0,10);
histos["nJets"] 	= new TH1F("nJets","nJets;N_{jets};events",10,0,10);

//rho
histos["rho"]		= new TH1F("rho","rho;#rho;events",500,0,50);
histos["rho25"]		= new TH1F("rho25","rho25;#rho(2.5); events",500,0,50);

//DiBoson
histos["llM"]		= new TH1F("llM","llM;M^{ll};events",500,0,150);
histos["llPt"]		= new TH1F("llPt","llPt;P_{T}^{ll};events",500,0,150);
histos["llY"]		= new TH1F("llY","llY;Y^{ll};events",500,-5,5);
histos["llPhi"]		= new TH1F("llPhi","llPhi;#phi^{ll};events",500,-3.1416,3.1416);

//photons?

//lepton
histos["lep0Pt"]	= new TH1F("lep0Pt","lep0Pt;P_{T}^{1st lep};events",500,0,150);
histos["lep0Eta"]	= new TH1F("lep0Eta","lep0Eta;#eta^{1st lep};events",500,-5,5);
histos["lep0Phi"]	= new TH1F("lep0Phi","lep0Phi;#phi^{1st lep};events",500,-3.1416,3.1416);

histos["lep1Pt"]	= new TH1F("lep1Pt","lep1Pt;#P_{T}^{2nd lep};events",500,0,150);
histos["lep1Eta"]	= new TH1F("lep1Eta","lep1Eta;#eta^{2nd lep};events",500,-5,5);
histos["lep1Phi"]	= new TH1F("lep1Phi","lep1Phi;#phi^{2nd lep};events",500,-3.1416,3.1416);

//jets
histos["jet0Pt"]	= new TH1F("jet0Pt","jet0Pt;P_{T}^{1st jet};events",500,0,150);
histos["jet0Eta"]	= new TH1F("jet0Eta","jet0Eta;#eta^{1st jet};events",500,-5,5);
histos["jet0Phi"]	= new TH1F("jet0Phi","jet0Phi;#phi^{1st jet};events",500,-3.1416,3.1416);
histos["jet0QGL"]	= new TH1F("jet0QGL","jet0QGL;QGL^{1st jet};events",500,-1.001,1.001);
histos["jet0Btag"]	= new TH1F("jet0Btag","jet0Btag;Btag^{1st jet};events",500,-1.001,1.001);

histos["jet1Pt"]	= new TH1F("jet1Pt","jet1Pt;P_{T}^{2nd jet};events",500,0,150);
histos["jet1Eta"]	= new TH1F("jet1Eta","jet1Eta;#eta^{2nd jet};events",500,-5,5);
histos["jet1Phi"]	= new TH1F("jet1Phi","jet1Phi;#phi^{2nd jet};events",500,-3.1416,3.1416);
histos["jet1QGL"]	= new TH1F("jet1QGL","jet1QGL;QGL^{1st jet};events",500,-1.001,1.001);
histos["jet1Btag"]	= new TH1F("jet1Btag","jet1Btag;Btag^{2nd jet};events",500,-1.001,1.001);
return 0;
}
//------------------- LOOP
void Loop(map<string,TH1F*> &histos, TTree *t,int type)
{
//DEBUG LEVEL 2
int selRECO;		t->SetBranchAddress("selRECO",&selRECO);
int nVtx;		t->SetBranchAddress("nVtx",&nVtx);
int nLeptons;		t->SetBranchAddress("nLeptons",&nLeptons);
int nPhotons;		t->SetBranchAddress("nPhotons",&nPhotons);
int nJets;		t->SetBranchAddress("nJets",&nJets);
float rho;		t->SetBranchAddress("rho",&rho);
float rho25;		t->SetBranchAddress("rho25",&rho25);
float llM;		t->SetBranchAddress("llM",&llM);
float llPt;		t->SetBranchAddress("llPt",&llPt);
float llPhi;		t->SetBranchAddress("llPhi",&llPhi);
float llY;		t->SetBranchAddress("llY",&llY);
float llDPhi;		t->SetBranchAddress("llDPhi",&llDPhi);
float llEta;		t->SetBranchAddress("llEta",&llEta);
vector<float> *lepPt=NULL	;t->SetBranchAddress("lepPt",&lepPt);
vector<float> *lepEta=NULL	;t->SetBranchAddress("lepEta",&lepEta);
vector<float> *lepPhi=NULL	;t->SetBranchAddress("lepPhi",&lepPhi);
vector<float> *lepE=NULL	;t->SetBranchAddress("lepE",&lepE);
vector<int>   *lepChId=NULL	;t->SetBranchAddress("lepChId",&lepChId);
vector<float> *jetPt=NULL	;t->SetBranchAddress("jetPt",&jetPt);
vector<float> *jetEta=NULL	;t->SetBranchAddress("jetEta",&jetEta);
vector<float> *jetPhi=NULL	;t->SetBranchAddress("jetPhi",&jetPhi);
vector<float> *jetE=NULL	;t->SetBranchAddress("jetE",&jetE);
vector<float> *jetQGL=NULL	;t->SetBranchAddress("jetQGL",&jetQGL);
vector<float> *jetBtag=NULL	;t->SetBranchAddress("jetBtag",&jetBtag);
vector<int>   *jetVeto=NULL	;t->SetBranchAddress("jetVeto",&jetVeto);

if(debug>1)cout<<"Beginning loop "<<endl;
for(unsigned long long iEntry=0;iEntry<t->GetEntries();iEntry++)
	{
	t->GetEntry(iEntry);
	if(debug>1)cout<<"Getting Entry "<<iEntry<<endl;
	if(!selRECO)continue;
	if(lepPt->size()<nLeptons)cout<<"WARINING nLeps"<<endl;
	if(lepPt->size()<2)continue; //2 leptons
	if( (*lepChId)[0]*(*lepChId)[1]!=-4)continue; //two muons
	//llM selection
	if(!(fabs(llM-91)<20))continue;
	//Triggers??
	if(debug>1)cout<<"Passed Selection"<<endl;
	int nJetsVeto=0;
	for(int iJ=0;iJ<jetVeto->size();iJ++)if( !( ( (*jetVeto)[iJ]&1)  && ( (*jetVeto)[iJ]&2)) ) nJetsVeto++;
	if(debug>1)cout<<"nJetsVeto "<<nJetsVeto<<endl;
	
	histos["nVtx"]->Fill(nVtx);
	histos["nLeptons"]->Fill(nLeptons);
	histos["nJets"]->Fill(nJetsVeto);

	histos["llM"]->Fill(llM);
	histos["llPt"]->Fill(llPt);
	histos["llY"]->Fill(llY);
	histos["llPhi"]->Fill(llPhi);
	
	if(debug>1)cout<<"lepSize Pt "<<lepPt->size()<<" Eta "<<lepEta->size()<<" Phi "<<lepPhi->size()<<endl;
	histos["lep0Pt"]	->Fill((*lepPt)[0]);
	histos["lep1Pt"]	->Fill((*lepPt)[1]);
	histos["lep0Eta"]	->Fill((*lepEta)[0]);
	histos["lep1Eta"]	->Fill((*lepEta)[1]);
	histos["lep0Phi"]	->Fill((*lepPhi)[0]);
	histos["lep1Phi"]	->Fill((*lepPhi)[1]);

	if(debug>1)cout<<"jetSize Pt "<<jetPt->size()<<" Veto "<<jetVeto->size()<<endl;
	int jet0=-1,jet1=-1;
		if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  && ( (*jetVeto)[iJ]&2)) ){jet0=iJ; break;}
		if(nJetsVeto>1)for(int iJ=jet0+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  && ((*jetVeto)[iJ]&2)) ){jet1=iJ; break;}
	if(debug>1)cout<<"Jet 0 Idx "<<jet0<<"Jet 1 Idx "<<jet1<<endl;

	if(jet0>=0)histos["jet0Pt"]	->Fill((*jetPt)[jet0]);
	if(jet1>=0)histos["jet1Pt"]	->Fill((*jetPt)[jet1]);
	if(jet0>=0)histos["jet0Eta"]	->Fill((*jetEta)[jet0]);
	if(jet1>=0)histos["jet1Eta"]	->Fill((*jetEta)[jet1]);
	if(jet0>=0)histos["jet0Phi"]	->Fill((*jetPhi)[jet0]);
	if(jet1>=0)histos["jet1Phi"]	->Fill((*jetPhi)[jet1]);
	if(jet0>=0)histos["jet0QGL"]	->Fill((*jetQGL)[jet0]);
	if(jet1>=0)histos["jet1QGL"]	->Fill((*jetQGL)[jet1]);
	if(jet0>=0)histos["jet0Btag"]	->Fill((*jetBtag)[jet0]);
	if(jet1>=0)histos["jet1Btag"]	->Fill((*jetBtag)[jet1]);
	}
}

//------------------- MAIN PROCEDURE ---------------------------------
void ValidationPlot(){
//DEBUG LEVEL 1
ValidationPlotData("/tmp/amarini/DoubleMu_Run2012A-13Jul2012-v1_AOD.root","Plot/validationHistosData.root","Plot/");
return; 
}

//------------------- Produce Plots DATA ---------------------------------
void ValidationPlotData(const char *fileIn,const char *fileOut, const char *dirOut){
//DEBUG LEVEL 1
map<string,TH1F*> histos;

if(debug>0)cout<<"Validation Plot Data"<<endl;
TFile *fOut=TFile::Open(fileOut,"UPDATE");
TFile *fIn=TFile::Open(fileIn);
TTree *tIn=(TTree*)fIn->Get("accepted/events");

if(debug>0)cout<<"Create Histos"<<endl;
CreateHisto(histos);
if(debug>0)cout<<"Begin Loop"<<endl;
Loop(histos,tIn,0);

if(debug>0)cout<<"Saving results"<<endl;
for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++){
	TCanvas *c=new TCanvas("c","c",800,800);
	iM->second->SetMarkerStyle(24);
	iM->second->Draw("P");
	c->SaveAs(Form("%s/%s",dirOut,iM->first.c_str()));
	}

fOut->cd();
for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++)
	iM->second->Write();
	
}

//------------------- MAIN -------------------------------------------
#ifdef STANDALONE
int main(){
//DEBUG LEVEL 1
}
#endif
