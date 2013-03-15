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
histos["llgM"]		= new TH1F("llgM","llgM;M^{ll#gamma};events/10GeV",200,0,2000);
histos["l1gM"]		= new TH1F("l1gM","l1gM;M^{l1#gamma};events/10GeV",200,0,2000);
histos["l2gM"]		= new TH1F("l2gM","l2gM;M^{l2#gamma};events/10GeV",200,0,2000);

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

// a bit of variables
histos["jetLLDPhi0"]	= new TH1F("JetLLDPhi0","JetLLDPhi0;#Delta#phi(Z,j_{1});events",50,0,3.1416);
histos["Sum3j"]	= new TH1F("Sum3j","Sum3j;#sum_{ij#elem 1..3}d#phi_{ij};events",50,0,3.1416*2);

//Introduce some cuts:
histos["llPt_betaStar"]		= new TH1F("llPt_betaStar","llPt;P_{T}^{ll} (betaStar on Jets);events",100,0,450);
histos["jetLLDPhi0_PtZ_50"]	= new TH1F("JetLLDPhi0_PtZ_50","JetLLDPhi0;#Delta#phi(Z,j_{1}) [P_{T}^{ll}>50 GeV];events",50,0,3.1416);
histos["llPt_nJets_3"]		= new TH1F("llPt_nJets_3","llPt;P_{T}^{ll} (N_{jets} #geq 3);events",100,0,450);

for(map<string,TH1F*>::iterator it=histos.begin();it!=histos.end();it++) it->second->Sumw2(); 

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
vector<float> *jetBeta=NULL	;t->SetBranchAddress("jetBeta",&jetBeta);
vector<int>   *jetVeto=NULL	;//t->SetBranchAddress("jetVeto",&jetVeto);
	jetVeto=new vector<int>;
vector<float> *photonPt=NULL; t->SetBranchAddress("photonPt",&photonPt);
vector<float> *photonEta=NULL; t->SetBranchAddress("photonEta",&photonEta);
vector<float> *photonPhi=NULL; t->SetBranchAddress("photonPhi",&photonPhi);
vector<float> *photonE=NULL; t->SetBranchAddress("photonE",&photonE);

//vector<float> *photonIsoFPRNeutral=NULL;if(type==0) t->SetBranchAddress("photonIsoFPRNeutral",&photonIsoFPRNeutral);
//vector<float> *photonIsoFPRCharged=NULL;if(type==0) t->SetBranchAddress("photonIsoFPRCharged",&photonIsoFPRCharged);
//vector<float> *photonIsoFPRPhoton=NULL; if(type==0) t->SetBranchAddress("photonIsoFPRPhoton",&photonIsoFPRPhoton);


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
		(*jetVeto)[k]=0;
		if(R1<JetDRCut) (*jetVeto)[k]=1;
		if(R2<JetDRCut) (*jetVeto)[k]|=2;

		}
		
	int nJetsVeto=0;
	for(int iJ=0;iJ<jetVeto->size() && (*jetPt)[iJ]>=JetPtCut;iJ++)if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2)) ) nJetsVeto++;
	if(debug>1)cout<<"nJetsVeto "<<nJetsVeto<<endl;
	if(debug>1)cout<<"jetSize Pt "<<jetPt->size()<<" Veto "<<jetVeto->size()<<endl;
	int jet0=-1,jet1=-1,jet2=-1;
		if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2)) ){jet0=iJ; break;}
		if(nJetsVeto>1)for(int iJ=jet0+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)) ){jet1=iJ; break;}
		if(nJetsVeto>2)for(int iJ=jet1+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)) ){jet2=iJ; break;}
	if( (jet0>=0) && (*jetPt)[jet0]< JetPtCut) jet0=-1;
	if( (jet1>=0) && (*jetPt)[jet1]< JetPtCut) jet1=-1;
	if( (jet2>=0) && (*jetPt)[jet2]< JetPtCut) jet2=-1;

	double weight=1;
	if(type>0)weight=PUWeight;


	if( jet0 < 0 )continue; //cut on the first jet

	if((photonPt->size()>0)&&((*photonPt)[0]>130 )){
		//find the first isolated photon
		int pho0=0;
		//if(type==0){for(pho0=0;pho0<photonPt->size();pho0++) if( ((*photonIsoFPRNeutral)[pho0]+(*photonIsoFPRCharged)[pho0]+(*photonIsoFPRPhoton)[pho0])/(*photonPt)[pho0] < .5) {break;} }
		TLorentzVector g;
		if( (pho0<photonPt->size()) && ((*photonPt)[pho0]>130)){
		g.SetPtEtaPhiE( (*photonPt)[pho0],(*photonEta)[pho0],(*photonPhi)[pho0],(*photonE)[pho0]);
		histos["llgM"]->Fill((l1+l2+g).M(),weight);
		histos["l1gM"]->Fill((l1+g).M(),weight);
		histos["l2gM"]->Fill((l2+g).M(),weight);
		}
	}
	//llM selection
	if(!(fabs(llM-91)<llMCut))continue;
	
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
	
	
	TLorentzVector ll=l1+l2,j0,j1,j2;
	bool isZlead=false;
		if( (jet1>=0) &&(llPt>(*jetPt)[jet1]))isZlead=true;
	if(jet0>=0) j0.SetPtEtaPhiE((*jetPt)[jet0],(*jetEta)[jet0],(*jetPhi)[jet0],(*jetE)[jet0]);
	if(jet1>=0) j1.SetPtEtaPhiE((*jetPt)[jet1],(*jetEta)[jet1],(*jetPhi)[jet1],(*jetE)[jet1]);
	if(jet2>=0) j2.SetPtEtaPhiE((*jetPt)[jet2],(*jetEta)[jet2],(*jetPhi)[jet2],(*jetE)[jet2]);

	if(jet0>=0){histos["jetLLDPhi0"]->Fill(fabs(j0.DeltaPhi(ll)),weight);}
	if(jet2>=0 && isZlead){histos["Sum3j"]->Fill( fabs(j0.DeltaPhi(j1))+fabs(j1.DeltaPhi(j2))+fabs(j0.DeltaPhi(j2)),weight);}

	//----------------------Jets with BStar: requirements already applied  E jet0 + veto
	for(int k=0;k<jetPt->size();k++)
		{
		if( (*jetPt)[k] <JetPtCut) continue;
		if( ( 1. - (*jetBeta)[k] >= 0.2 * TMath::Log( nVtx - 0.67))) //betaStar=1-beta
		 	(*jetVeto)[k]|=8; //1= lept1, 2=lept2, 4= gamma, 8=betaStar
		}
	 nJetsVeto=0;
	for(int iJ=0;iJ<jetVeto->size() && (*jetPt)[iJ]>=JetPtCut;iJ++)if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2)) ) nJetsVeto++;

	int jet0_BS=-1,jet1_BS=-1,jet2_BS=-1;
		if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2) || ((*jetVeto)[iJ]&8)) ) {jet0_BS=iJ; break;}
		if(nJetsVeto>1)for(int iJ=jet0_BS+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)|| ((*jetVeto)[iJ]&8) ) ){jet1_BS=iJ; break;}
		if(nJetsVeto>2)for(int iJ=jet1_BS+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)|| ((*jetVeto)[iJ]&8)) ){jet2_BS=iJ; break;}

	if( (jet0_BS>=0) && (*jetPt)[jet0_BS]< JetPtCut) jet0_BS=-1;
	if( (jet1_BS>=0) && (*jetPt)[jet1_BS]< JetPtCut) jet1_BS=-1;
	if( (jet2_BS>=0) && (*jetPt)[jet2_BS]< JetPtCut) jet2_BS=-1;

	if(jet0_BS>=0)histos["llPt_betaStar"]->Fill(llPt,weight);
        if( (jet0>=0) && (llPt>50) )   histos["jetLLDPhi0_PtZ_50"]->Fill(fabs(j0.DeltaPhi(ll)),weight);
	if( (jet2>=0) )histos["llPt_nJets_3"]->Fill(llPt,weight);	

	} //END LOOP OVER ENTRIES
	
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
//TFile *fIn=TFile::Open(fileIn); //TTree *t=(TTree*)fIn->Get("accepted/events");
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
sscanf(A.ReadParameter("BATCH"),"%d",&useFork)  ;
//not in the configuration fail
int forkNum=0;
if(useFork)
	{
	if(argc<3)useFork=0;
	else sscanf(argv[2],"%d",&forkNum);
	}
	

if( useFork){

	if(forkNum==1)ValidationPlotMC( (DirMC+DY).c_str(),Form("DY_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==2)ValidationPlotMC( (DirMC+TT).c_str(),Form("TT_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==3)ValidationPlotMC( (DirMC+WJ).c_str(),Form("WJ_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==4)ValidationPlotMC( (DirMC+WW).c_str(),Form("WW_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==5)ValidationPlotMC( (DirMC+WZ).c_str(),Form("WZ_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==6)ValidationPlotMC( (DirMC+ZZ).c_str(),Form("ZZ_%d.root",-CHID2), DirOut.c_str());
	
	if(CHID2==-4){if(forkNum==7)ValidationPlotData( (DirData+DoubleMu).c_str(),"DoubleMu_4.root", DirOut.c_str());}
	if(CHID2==-1){if(forkNum==7)ValidationPlotData( (DirData+DoubleE).c_str(),"DoubleE_1.root", DirOut.c_str());}
	if(CHID2==-2){if(forkNum==7)ValidationPlotData( (DirData+MuEG).c_str(),"MuEG_2.root", DirOut.c_str());}

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
