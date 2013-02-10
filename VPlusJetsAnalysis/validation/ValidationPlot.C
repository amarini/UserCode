#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <unistd.h>

#ifdef STANDALONE
	#include "src/ReadParameters.C"
#endif

using namespace std;
const int debug=1;//0 NO DEBUG 1 MINIMAL 2 MAX

void ValidationPlotData(const char *fileIn,const char *fileOut, const char *dirOut);
int CreateHisto(map<string,TH1F*> &histos);
void ValidationPlot();
void Loop(map<string,TH1F*> &histos, TChain *t,int type=0); //0 data 1 mc

float JetPtCut=50;
float JetDRCut=0.4;
float llMCut=20;
int CHID2=-4;


// ------------------ Varibales plots -------------------------
int CreateHisto(map<string,TH1F*> &histos)
{
for(map<string,TH1F*>::iterator it=histos.begin();it!=histos.end();it++)
	{
	delete it->second;
	}
histos.clear();
//DEBUG LEVEL 2
//--- counting
histos["nVtx"] 		= new TH1F("nVtx","nVtx;N_{Vtx};events",50,0,50);
histos["nLeptons"] 	= new TH1F("nLeptons","nLeptons;N_{lep};events",10,0,10);
histos["nPhotons"] 	= new TH1F("nPhotons","nPhotons;N_{#gamma};events",10,0,10);
histos["nJets"] 	= new TH1F("nJets","nJets;N_{jets};events",10,0,10);

//rho
histos["rho"]		= new TH1F("rho","rho;#rho;events",50,0,50);
histos["rho25"]		= new TH1F("rho25","rho25;#rho(2.5); events",50,0,50);

//DiBoson
histos["llM"]		= new TH1F("llM","llM;M^{ll};events",50,40,150);
histos["llPt"]		= new TH1F("llPt","llPt;P_{T}^{ll};events",100,0,450);
histos["llY"]		= new TH1F("llY","llY;Y^{ll};events",50,-5,5);
histos["llPhi"]		= new TH1F("llPhi","llPhi;#phi^{ll};events",50,-3.1416,3.1416);

//photons?

//lepton
histos["lep0Pt"]	= new TH1F("lep0Pt","lep0Pt;P_{T}^{1st lep};events",50,20,150);
histos["lep0Eta"]	= new TH1F("lep0Eta","lep0Eta;#eta^{1st lep};events",50,-5,5);
histos["lep0Phi"]	= new TH1F("lep0Phi","lep0Phi;#phi^{1st lep};events",50,-3.1416,3.1416);

histos["lep1Pt"]	= new TH1F("lep1Pt","lep1Pt;#P_{T}^{2nd lep};events",50,20,150);
histos["lep1Eta"]	= new TH1F("lep1Eta","lep1Eta;#eta^{2nd lep};events",50,-5,5);
histos["lep1Phi"]	= new TH1F("lep1Phi","lep1Phi;#phi^{2nd lep};events",50,-3.1416,3.1416);

//jets
histos["jet0Pt"]	= new TH1F("jet0Pt","jet0Pt;P_{T}^{1st jet};events",50,50,150);
histos["jet0Eta"]	= new TH1F("jet0Eta","jet0Eta;#eta^{1st jet};events",50,-5,5);
histos["jet0Phi"]	= new TH1F("jet0Phi","jet0Phi;#phi^{1st jet};events",50,-3.1416,3.1416);
histos["jet0QGL"]	= new TH1F("jet0QGL","jet0QGL;QGL^{1st jet};events",50,-1.001,1.001);
histos["jet0Btag"]	= new TH1F("jet0Btag","jet0Btag;Btag^{1st jet};events",50,-1.001,1.001);

histos["jet1Pt"]	= new TH1F("jet1Pt","jet1Pt;P_{T}^{2nd jet};events",50,50,150);
histos["jet1Eta"]	= new TH1F("jet1Eta","jet1Eta;#eta^{2nd jet};events",50,-5,5);
histos["jet1Phi"]	= new TH1F("jet1Phi","jet1Phi;#phi^{2nd jet};events",50,-3.1416,3.1416);
histos["jet1QGL"]	= new TH1F("jet1QGL","jet1QGL;QGL^{1st jet};events",50,-1.001,1.001);
histos["jet1Btag"]	= new TH1F("jet1Btag","jet1Btag;Btag^{2nd jet};events",50,-1.001,1.001);
return 0;
}
//------------------- LOOP
void Loop(map<string,TH1F*> &histos, TChain *t,int type)
{
//DEBUG LEVEL 2
//int selRECO;		t->SetBranchAddress("selRECO",&selRECO);
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
vector<int>   *jetVeto=NULL	;//t->SetBranchAddress("jetVeto",&jetVeto);
	jetVeto=new vector<int>;

double PUWeight; if(type>0)t->SetBranchAddress("PUWeight",&PUWeight);


if(debug>1)cout<<"Beginning loop "<<endl;
for(unsigned long long iEntry=0;iEntry<t->GetEntries() ;iEntry++)
	{
	t->GetEntry(iEntry);
	if(debug>1)cout<<"Getting Entry "<<iEntry<<endl;
	//if(!selRECO)continue;
	//if(lepPt->size()<nLeptons)cout<<"WARNING nLeps"<<endl;
	if(lepPt->size()<2)continue; //2 leptons
	if( (*lepChId)[0]*(*lepChId)[1]!=CHID2)continue; //two muons
	//llM selection
	if(!(fabs(llM-91)<llMCut))continue;
	//Triggers??
	if(debug>1)cout<<"Passed Selection"<<endl;
		//redo vetos
	jetVeto->clear();
	jetVeto->resize(jetPt->size());
	TLorentzVector j,l1,l2;
	l1.SetPtEtaPhiE( (*lepPt)[0],(*lepEta)[0],(*lepPhi)[0],(*lepE)[0]);
	l2.SetPtEtaPhiE( (*lepPt)[1],(*lepEta)[1],(*lepPhi)[1],(*lepE)[1]);
	for(int k=0;k<jetPt->size();k++)
		{
		if( (*jetPt)[k] <JetPtCut) continue;
			j.SetPtEtaPhiE( (*jetPt)[k],(*jetEta)[k],(*jetPhi)[k],(*jetE)[k]);
		float R1=l1.DeltaR(j);
		float R2=l2.DeltaR(j);
		if(R1<JetDRCut) (*jetVeto)[k]=1;
		if(R2<JetDRCut) (*jetVeto)[k]|=2;
		else (*jetVeto)[k]=0;

		}
		
	int nJetsVeto=0;
	for(int iJ=0;iJ<jetVeto->size() && (*jetPt)[iJ]>=JetPtCut;iJ++)if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2)) ) nJetsVeto++;
	if(debug>1)cout<<"nJetsVeto "<<nJetsVeto<<endl;

	double weight=1;
	if(type>0)weight=PUWeight;
	
	histos["nVtx"]->Fill(nVtx,weight);
	histos["nLeptons"]->Fill(nLeptons,weight);
	histos["nPhotons"]->Fill(nPhotons,weight);
	histos["nJets"]->Fill(nJetsVeto,weight);

	histos["rho"]->Fill(rho,weight);
	histos["rho25"]->Fill(rho25,weight);

	histos["llM"]->Fill(llM,weight);
	histos["llPt"]->Fill(llPt,weight);
	histos["llY"]->Fill(llY,weight);
	histos["llPhi"]->Fill(llPhi,weight);
	
	if(debug>1)cout<<"lepSize Pt "<<lepPt->size()<<" Eta "<<lepEta->size()<<" Phi "<<lepPhi->size()<<endl;
	histos["lep0Pt"]	->Fill((*lepPt)[0],weight);
	histos["lep1Pt"]	->Fill((*lepPt)[1],weight);
	histos["lep0Eta"]	->Fill((*lepEta)[0],weight);
	histos["lep1Eta"]	->Fill((*lepEta)[1],weight);
	histos["lep0Phi"]	->Fill((*lepPhi)[0],weight);
	histos["lep1Phi"]	->Fill((*lepPhi)[1],weight);

	if(debug>1)cout<<"jetSize Pt "<<jetPt->size()<<" Veto "<<jetVeto->size()<<endl;
	int jet0=-1,jet1=-1;
		if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2)) ){jet0=iJ; break;}
		if(nJetsVeto>1)for(int iJ=jet0+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)) ){jet1=iJ; break;}
	if(debug>1)cout<<"Jet 0 Idx "<<jet0<<"Jet 1 Idx "<<jet1<<endl;

	if(jet0>=0)histos["jet0Pt"]	->Fill((*jetPt)[jet0],weight);
	if(jet1>=0)histos["jet1Pt"]	->Fill((*jetPt)[jet1],weight);
	if(jet0>=0)histos["jet0Eta"]	->Fill((*jetEta)[jet0],weight);
	if(jet1>=0)histos["jet1Eta"]	->Fill((*jetEta)[jet1],weight);
	if(jet0>=0)histos["jet0Phi"]	->Fill((*jetPhi)[jet0],weight);
	if(jet1>=0)histos["jet1Phi"]	->Fill((*jetPhi)[jet1],weight);
	if(jet0>=0)histos["jet0QGL"]	->Fill((*jetQGL)[jet0],weight);
	if(jet1>=0)histos["jet1QGL"]	->Fill((*jetQGL)[jet1],weight);
	if(jet0>=0)histos["jet0Btag"]	->Fill((*jetBtag)[jet0],weight);
	if(jet1>=0)histos["jet1Btag"]	->Fill((*jetBtag)[jet1],weight);
	}
	jetVeto->clear();
	delete jetVeto;
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
//TFile *fOut=TFile::Open(fileOut,"UPDATE");
TFile *fOut=TFile::Open(  Form("%s%s",dirOut,fileOut),"UPDATE");
//TFile *fIn=TFile::Open(fileIn);
TChain *tIn=new TChain("accepted/events");
Int_t nFiles=tIn->Add(fileIn);
cout<<"Added "<<nFiles <<" files to the chain"<<endl;
if(debug>0) for(int i=0;i<nFiles;i++) cout<<"   -->File "<<i<<" is " <<tIn->GetListOfFiles()->At(i)->GetTitle()<<endl;
//TTree *tIn=(TTree*)fIn->Get("accepted/events");

if(debug>0)cout<<"Create Histos"<<endl;
CreateHisto(histos);
if(debug>0)cout<<"Begin Loop"<<endl;
Loop(histos,tIn,0);

if(debug>0)cout<<"Saving results"<<endl;
//for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++){
//	TCanvas *c=new TCanvas("c","c",800,800);
//	iM->second->SetMarkerStyle(24);
//	iM->second->Draw("P");
//	c->SaveAs(Form("%s/%s.pdf",dirOut,iM->first.c_str()));
//	}

fOut->cd();
for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++)
	iM->second->Write();

fOut->Close();
//for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++)iM->second->Delete();
histos.clear();
tIn->Delete();
	
}

//------------------- Produce Plots MC ---------------------------------
void ValidationPlotMC(const char *fileIn,const char *fileOut, const char *dirOut){
//DEBUG LEVEL 1
map<string,TH1F*> histos;

if(debug>0)cout<<"Validation Plot MC"<<endl;
TFile *fOut=TFile::Open(  Form("%s%s",dirOut,fileOut),"RECREATE");
TFile *fIn=TFile::Open(fileIn); //TTree *t=(TTree*)fIn->Get("accepted/events");
TChain *tIn=new TChain("accepted/events");
Int_t nFiles=tIn->Add(fileIn);
cout<<"Added "<< nFiles <<" files to the chain"<<endl;
if(debug>0) for(int i=0;i<nFiles;i++) cout<<"   -->File "<<i<<" is " <<tIn->GetListOfFiles()->At(i)->GetTitle()<<endl;

if(debug>0)cout<<"Create Histos"<<endl;
CreateHisto(histos);
if(debug>0)cout<<"Begin Loop"<<endl;
Loop(histos,tIn,1);

if(debug>0)cout<<"Saving results"<<endl;
//for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++){
//	TCanvas *c=new TCanvas("c","c",800,800);
//	iM->second->SetMarkerStyle(24);
//	iM->second->Draw("P");
//	c->SaveAs(Form("%s/%s.pdf",dirOut,iM->first.c_str()));
//	}

fOut->cd();
for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++)
	iM->second->Write();
fOut->Close();
//for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++)iM->second->Delete();
histos.clear();
tIn->Delete();
}

//------------------- MAIN -------------------------------------------
#ifdef STANDALONE
int main(int argc, char **argv){
//DEBUG LEVEL 1
string configFile;

if(argc<2) configFile="data/config.ini";
else configFile=argv[1];

Read A(configFile.c_str());

sscanf(A.ReadParameter("JETPT"),"%f",&JetPtCut );
sscanf(A.ReadParameter("JETDR"),"%f",&JetDRCut);
sscanf(A.ReadParameter("LLM"),"%f",&llMCut );
sscanf(A.ReadParameter("CHID2"),"%d",&CHID2 );

printf("********CUT********\n");
printf("* JetPt=%4.1f      *\n",JetPtCut);
printf("* JetDR=%4.2f      *\n",JetDRCut);
printf("* llM  =%4.1f      *\n",llMCut);
printf("* ChID^2=%3d      *\n",CHID2);
printf("*******************\n");

string DirMC=A.ReadParameter("MCDIR"); DirMC+="/";
string DY=A.ReadParameter("DY");
string TT=A.ReadParameter("TT");
string WJ=A.ReadParameter("WJ");
string WW=A.ReadParameter("WW");
string WZ=A.ReadParameter("WZ");
string ZZ=A.ReadParameter("ZZ");

string DirData=A.ReadParameter("DATADIR"); DirData+="/";
string DoubleMu=A.ReadParameter("DoubleMu");
string DoubleE=A.ReadParameter("DoubleE");
string MuEG=A.ReadParameter("MuEG");
 
string DirOut=A.ReadParameter("OUTDIR"); DirOut+="/";

int useFork=0;
if( useFork){

	int who=fork(); //pid son, 0 parent
	if(who==0)ValidationPlotMC( (DirMC+DY).c_str(),Form("DY_%d.root",-CHID2), DirOut.c_str());
	if(who==0)return 0;
	who=fork();if(who==0)ValidationPlotMC( (DirMC+TT).c_str(),Form("TT_%d.root",-CHID2), DirOut.c_str());
	if(who==0)return 0;
	who=fork();if(who==0)ValidationPlotMC( (DirMC+WJ).c_str(),Form("WJ_%d.root",-CHID2), DirOut.c_str());
	if(who==0)return 0;
	who=fork();if(who==0)ValidationPlotMC( (DirMC+WW).c_str(),Form("WW_%d.root",-CHID2), DirOut.c_str());
	if(who==0)return 0;
	who=fork();if(who==0)ValidationPlotMC( (DirMC+WZ).c_str(),Form("WZ_%d.root",-CHID2), DirOut.c_str());
	if(who==0)return 0;
	who=fork();if(who==0)ValidationPlotMC( (DirMC+ZZ).c_str(),Form("ZZ_%d.root",-CHID2), DirOut.c_str());
	if(who==0)return 0;
	
	if(CHID2==-4){who=fork();if(who==0)ValidationPlotData( (DirData+DoubleMu).c_str(),"DoubleMu", DirOut.c_str());}
	if(who==0)return 0;
	if(CHID2==-1){who=fork();if(who==0)ValidationPlotData( (DirData+DoubleE).c_str(),"DoubleE", DirOut.c_str());}
	if(who==0)return 0;
	if(CHID2==-2){who=fork();if(who==0)ValidationPlotData( (DirData+MuEG).c_str(),"MuEG", DirOut.c_str());}
	if(who==0)return 0;

} else { //!FORK
	ValidationPlotMC( (DirMC+DY).c_str(),Form("DY_%d.root",-CHID2), DirOut.c_str());
	ValidationPlotMC( (DirMC+TT).c_str(),Form("TT_%d.root",-CHID2), DirOut.c_str());
	ValidationPlotMC( (DirMC+WJ).c_str(),Form("WJ_%d.root",-CHID2), DirOut.c_str());
	ValidationPlotMC( (DirMC+WW).c_str(),Form("WW_%d.root",-CHID2), DirOut.c_str());
	ValidationPlotMC( (DirMC+WZ).c_str(),Form("WZ_%d.root",-CHID2), DirOut.c_str());
	ValidationPlotMC( (DirMC+ZZ).c_str(),Form("ZZ_%d.root",-CHID2), DirOut.c_str());
	
	if(CHID2==-4){ValidationPlotData( (DirData+DoubleMu).c_str(),"DoubleMu_4.root", DirOut.c_str());}
	if(CHID2==-1){ValidationPlotData( (DirData+DoubleE).c_str(),"DoubleE_1.root", DirOut.c_str());}
	if(CHID2==-2){ValidationPlotData( (DirData+MuEG).c_str(),"MuEG_2.root", DirOut.c_str());}

}
}
#endif
