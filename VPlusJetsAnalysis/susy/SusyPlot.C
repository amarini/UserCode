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

void SusyPlotData(const char *fileIn,const char *fileOut, const char *dirOut);
int CreateHisto(map<string,TH1F*> &histos);
void SusyPlot(){}
void Loop(map<string,TH1F*> &histos, TChain *t,int type=0); //0 data 1 mc

float JetPtCut=50;
float JetDRCut=0.4;
float llMCut=20;
float pfMetCut=100;
float LepEtaCut=1.4;
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

//DiBoson
histos["llM"]		= new TH1F("llM","llM;M^{ll};events",80,20,220);
histos["llM_2j"]	= new TH1F("llM_2j","llM;M^{ll};events",80,20,220);
histos["llM_beforeCut"]	= new TH1F("llM_beforeCut","llM;M^{ll};events",80,20,220);

histos["llM_BS"]		= new TH1F("llM_BS","llM;M^{ll};events",80,20,220);
histos["llM_2j_BS"]		= new TH1F("llM_2j_BS","llM;M^{ll};events",80,20,220);
histos["llM_beforeCut_BS"]	= new TH1F("llM_beforeCut_BS","llM;M^{ll};events",80,20,220);

//photons?

//lepton

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
float pfmet;		t->SetBranchAddress("pfmet",&pfmet);
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
	if(llM<20)continue;//minimum requirements

	if( fabs((*lepEta)[0]) >LepEtaCut) continue;
	if( fabs((*lepEta)[1]) >LepEtaCut) continue;
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
	//compute jets
	int jet0=-1,jet1=-1,jet2=-1;
		if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2)) ){jet0=iJ; break;}
		if(nJetsVeto>1)for(int iJ=jet0+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)) ){jet1=iJ; break;}
		if(nJetsVeto>2)for(int iJ=jet1+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)) ){jet2=iJ; break;}
	if( (jet0>=0) && (*jetPt)[jet0]< JetPtCut) jet0=-1;
	if( (jet1>=0) && (*jetPt)[jet1]< JetPtCut) jet1=-1;
	if( (jet2>=0) && (*jetPt)[jet2]< JetPtCut) jet2=-1;
	
	//compute jets BetaStar
	int jet0_BS=-1,jet1_BS=-1,jet2_BS=-1;
		if(nJetsVeto>0)for(int iJ=0;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ( (*jetVeto)[iJ]&2) || ((*jetVeto)[iJ]&8)) ) {jet0_BS=iJ; break;}
		if(nJetsVeto>1)for(int iJ=jet0_BS+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)|| ((*jetVeto)[iJ]&8) ) ){jet1_BS=iJ; break;}
		if(nJetsVeto>2)for(int iJ=jet1_BS+1;iJ<jetPt->size();iJ++) if( !( ( (*jetVeto)[iJ]&1)  || ((*jetVeto)[iJ]&2)|| ((*jetVeto)[iJ]&8)) ){jet2_BS=iJ; break;}

	if( (jet0_BS>=0) && (*jetPt)[jet0_BS]< JetPtCut) jet0_BS=-1;
	if( (jet1_BS>=0) && (*jetPt)[jet1_BS]< JetPtCut) jet1_BS=-1;
	if( (jet2_BS>=0) && (*jetPt)[jet2_BS]< JetPtCut) jet2_BS=-1;


	double weight=1;
	if(type>0)weight=PUWeight;

	if( jet0 < 0 )continue; //cut on the first jet

	histos["llM_beforeCut"]->Fill(llM,weight);

	if( pfmet >= pfMetCut ){
		if( jet1 >=0 ) {//2j
	
			histos["llM_2j"]->Fill(llM,weight);
			}// 2j
	
		if( jet2 >= 0 ){ //cut on the third jet
			histos["llM"]->Fill(llM,weight);
			}//3j
	} //end MET Cut
	//--------------------------betaStar jets----------------
	if( jet0_BS <0 ) continue;
	histos["llM_beforeCut_BS"]->Fill(llM,weight);

	if( pfmet >= pfMetCut ){
		if( jet1_BS >=0 ) {//2j
	
			histos["llM_2j_BS"]->Fill(llM,weight);
			}// 2j
	
		if( jet2_BS >= 0 ){ //cut on the third jet
			histos["llM_BS"]->Fill(llM,weight);
			}//3j
	} //end MET Cut
	
	
	TLorentzVector ll=l1+l2,j0,j1,j2;

	}
	jetVeto->clear();
	delete jetVeto;
}


//------------------- Produce Plots DATA ---------------------------------
void SusyPlotData(const char *fileIn,const char *fileOut, const char *dirOut){
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

fOut->cd();
for(map<string,TH1F*>::iterator iM=histos.begin();iM!=histos.end();iM++)
	iM->second->Write();

fOut->Close();
histos.clear();
tIn->Delete();
	
}

//------------------- Produce Plots MC ---------------------------------
void SusyPlotMC(const char *fileIn,const char *fileOut, const char *dirOut){
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

string configFile2;
if(argc>=3){configFile2=argv[2];}
else configFile2="";

Read A(configFile.c_str());

	sscanf(A.ReadParFromMultFile(configFile2.c_str(),"JETPT"),"%f",&JetPtCut );
	sscanf(A.ReadParFromMultFile(configFile2.c_str(),"JETDR"),"%f",&JetDRCut);
	//sscanf(A.ReadParFromMultFile(configFile2.c_str(),"LLM"),"%f",&llMCut ); //Met Cut is minimal
	sscanf(A.ReadParFromMultFile(configFile2.c_str(),"CHID2"),"%d",&CHID2 );
	sscanf(A.ReadParFromMultFile(configFile2.c_str(),"METCUT"),"%f",&pfMetCut );
	sscanf(A.ReadParFromMultFile(configFile2.c_str(),"LEPETA"),"%f",&LepEtaCut );


printf("********CUT********\n");
printf("* JetPt=%4.1f      *\n",JetPtCut);
printf("* JetDR=%4.2f      *\n",JetDRCut);
printf("* LepEta=%4.1f      *\n",LepEtaCut);
printf("* ChID^2=%3d      *\n",CHID2);
printf("* MET  =%4.1f      *\n",pfMetCut);
printf("*******************\n");

string DirMC=A.ReadParFromMultFile(configFile2.c_str(),"MCDIR"); DirMC+="/";
string DY=A.ReadParFromMultFile(configFile2.c_str(),"DY");
string TT=A.ReadParFromMultFile(configFile2.c_str(),"TT");
string WJ=A.ReadParFromMultFile(configFile2.c_str(),"WJ");
string WW=A.ReadParFromMultFile(configFile2.c_str(),"WW");
string WZ=A.ReadParFromMultFile(configFile2.c_str(),"WZ");
string ZZ=A.ReadParFromMultFile(configFile2.c_str(),"ZZ");

string DirData=A.ReadParFromMultFile(configFile2.c_str(),"DATADIR"); DirData+="/";
string DoubleMu=A.ReadParFromMultFile(configFile2.c_str(),"DoubleMu");
string DoubleE=A.ReadParFromMultFile(configFile2.c_str(),"DoubleE");
string MuEG=A.ReadParFromMultFile(configFile2.c_str(),"MuEG");
 
string DirOut=A.ReadParFromMultFile(configFile2.c_str(),"OUTDIR"); DirOut+="/";

int useFork=0;
sscanf(A.ReadParFromMultFile(configFile2.c_str(),"BATCH"),"%d",&useFork)  ;
int forkNum=0;
if(useFork)
	{
	if(argc<4)useFork=0;
	else sscanf(argv[3],"%d",&forkNum);
	}

if( useFork){
	if(forkNum==1)SusyPlotMC( (DirMC+DY).c_str(),Form("Susy_DY_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==2)SusyPlotMC( (DirMC+TT).c_str(),Form("Susy_TT_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==3)SusyPlotMC( (DirMC+WJ).c_str(),Form("Susy_WJ_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==4)SusyPlotMC( (DirMC+WW).c_str(),Form("Susy_WW_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==5)SusyPlotMC( (DirMC+WZ).c_str(),Form("Susy_WZ_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==6)SusyPlotMC( (DirMC+ZZ).c_str(),Form("Susy_ZZ_%d.root",-CHID2), DirOut.c_str());
	
	if(CHID2==-4){if(forkNum==7)SusyPlotData( (DirData+DoubleMu).c_str(),"Susy_DoubleMu_4.root", DirOut.c_str());}
	if(CHID2==-1){if(forkNum==7)SusyPlotData( (DirData+DoubleE).c_str(),"Susy_DoubleE_1.root", DirOut.c_str());}
	if(CHID2==-2){if(forkNum==7)SusyPlotData( (DirData+MuEG).c_str(),"Susy_MuEG_2.root", DirOut.c_str());}

} else { //!FORK
	SusyPlotMC( (DirMC+DY).c_str(),Form("Susy_DY_%d.root",-CHID2), DirOut.c_str());
	SusyPlotMC( (DirMC+TT).c_str(),Form("Susy_TT_%d.root",-CHID2), DirOut.c_str());
	SusyPlotMC( (DirMC+WJ).c_str(),Form("Susy_WJ_%d.root",-CHID2), DirOut.c_str());
	SusyPlotMC( (DirMC+WW).c_str(),Form("Susy_WW_%d.root",-CHID2), DirOut.c_str());
	SusyPlotMC( (DirMC+WZ).c_str(),Form("Susy_WZ_%d.root",-CHID2), DirOut.c_str());
	SusyPlotMC( (DirMC+ZZ).c_str(),Form("Susy_ZZ_%d.root",-CHID2), DirOut.c_str());
	
	if(CHID2==-4){SusyPlotData( (DirData+DoubleMu).c_str(),"Susy_DoubleMu_4.root", DirOut.c_str());}
	if(CHID2==-1){SusyPlotData( (DirData+DoubleE).c_str(),"Susy_DoubleE_1.root", DirOut.c_str());}
	if(CHID2==-2){SusyPlotData( (DirData+MuEG).c_str(),"Susy_MuEG_2.root", DirOut.c_str());}

}
}
#endif
