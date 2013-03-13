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

void QGPlotData(const char *fileIn,const char *fileOut, const char *dirOut);
int CreateHisto(map<string,TH1F*> &histos);
void QGPlot(){}
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
if(debug>1) cout<<"Starting Creating histos"<<endl;

//DiBoson
histos["llM"]		= new TH1F("llM","llM;M^{ll};events",80,20,220);
histos["llPt"]		= new TH1F("llPt","llPt;P_{T}^{ll};events",100,0,500);
for(int i=0;i<=histos["llPt"]->GetNbinsX()+1;i++)
	{
		//QGL
		{
		string name=Form("qgl_llPt_bin%d",i);
		histos[name.c_str()] =new TH1F(name.c_str(),"qgl;QGL;events",100,-1.0,1.000000001);
		for(int j=0;j<7;j++)	
			{
			string name=Form("qgl_llPt_bin%d_flavor%d",i,(j==6)?21:j);
			histos[name.c_str()] =new TH1F(name.c_str(),"qgl;QGL;events",100,-1.0,1.000000001);
			}
		}
		//BTAG
		{
		string name=Form("btag_llPt_bin%d",i);
		histos[name.c_str()] =new TH1F(name.c_str(),"btag;BTAG;events",100,-1.0,1.000000001);
		for(int j=0;j<7;j++)	
			{
			string name=Form("btag_llPt_bin%d_flavor%d",i,(j==6)?21:j);
			histos[name.c_str()] =new TH1F(name.c_str(),"btag;BTAG;events",100,-1.0,1.000000001);
			}
		}
	}
//Flavor PT
for(int j=0;j<7;j++)	
	{
	string name=Form("llPt_flavor%d",(j==6)?21:j);
	histos[name.c_str()] =new TH1F(name.c_str(),"llpt;llPt(F);events",100,0,500.0);
	}

//photons?

//lepton

for(map<string,TH1F*>::iterator it=histos.begin();it!=histos.end();it++) it->second->Sumw2(); 
if(debug>1) cout<<"Histos created"<<endl;
if(debug>1)for(map<string,TH1F*>::iterator it=histos.begin();it!=histos.end();it++) cout<<"Created histo"<<it->first<<endl;; 

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
vector<int>   *jetVeto=NULL	;//t->SetBranchAddress("jetVeto",&jetVeto);
	jetVeto=new vector<int>;

double PUWeight; if(type>0)t->SetBranchAddress("PUWeight",&PUWeight);
//V00-06
	//vector<float> *QGVars=NULL	;if(type>0)t->SetBranchAddress("QGVars",&QGVars); 
//V00-07
	vector<int> *jetPdgId=NULL	;if(type>0)t->SetBranchAddress("jetPdgId",&jetPdgId); 
	vector<int> *jetIdGEN=NULL	;if(type>0)t->SetBranchAddress("jetIdGEN",&jetIdGEN); 


if(debug>1)cout<<"Beginning loop "<<endl;
for(unsigned long long iEntry=0;iEntry<t->GetEntries() ;iEntry++)
	{
	if(debug>1 && iEntry <100) cout<<"Processing entry"<<iEntry << " of "<<t->GetEntries()<<endl;

	t->GetEntry(iEntry);
	if(debug>1)cout<<"Getting Entry "<<iEntry<<endl;
	//if(!selRECO)continue;
	//if(lepPt->size()<nLeptons)cout<<"WARNING nLeps"<<endl;
	if(lepPt->size()<2)continue; //2 leptons
	if( (*lepChId)[0]*(*lepChId)[1]!=CHID2)continue; //two muons
	//llM selection
	if(llM<20)continue;//minimum requirements
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
	
	if((debug>1) &&(jetQGL->size()!=jetPt->size()))cout<<"QGL vector size does not match"<<endl;
	if((debug>1) &&(jetBtag->size()!=jetPt->size()))cout<<"QGL vector size does not match"<<endl;

	if(debug>1)cout<<"I'm here!"<<endl;

	histos["llM"]->Fill(llM,weight);
	histos["llPt"] ->Fill(llPt,weight);
//for(int i=1;i<=histos["llPt"]->GetNbinsX();i++)
//	{
	int binllPt=histos["llPt"]->FindBin(llPt);
	{
	string name=Form("qgl_llPt_bin%d",binllPt);
	if(debug>1)cout<<"I'm here!QGL "<<name<<"trying to access jet0="<<jet0<<"/"<<jetQGL->size()<<"histo"<<histos[name.c_str()]<<endl;
	histos[name.c_str()] ->Fill(  (*jetQGL)[jet0],weight );
	if(debug>1)cout<<"I'm here!BTAG"<<endl;
	name=Form("btag_llPt_bin%d",binllPt);
	histos[name.c_str()] ->Fill(  (*jetBtag)[jet0],weight );
	}
	if(type>0) { //mc
	//V00-06
	//int Jet0Flavor=fabs((*QGVars)[10*jet0]);
	//V00-07
	if(debug>1)cout<<"I'm here!PDGID "<<jet0<<"/"<<jetPdgId->size()<<endl;
	int Jet0Flavor=fabs((*jetPdgId)[jet0]);

	if(debug>1)cout<<"pdgId= "<<jetPdgId->at(jet0)<<endl;
	//LLP - FLAVOR	
	{
	string name=Form("llPt_flavor%d",Jet0Flavor);
	histos[name.c_str()] ->Fill(llPt,weight);
	}

	//QGL
	if(debug>1)cout<<"I'm here QGL!"<<endl;
		{
		string name=Form("qgl_llPt_bin%d_flavor%d",binllPt,Jet0Flavor);
		histos[name.c_str()] ->Fill((*jetQGL)[jet0],weight);
		}
	//BTAG
	if(debug>1)cout<<"I'm here BTAG!"<<endl;
		{
		string name=Form("btag_llPt_bin%d_flavor%d",binllPt,Jet0Flavor) ; 
		histos[name.c_str()] ->Fill( (*jetBtag)[jet0],weight);
		}
	}
	
	if(debug>1)cout<<"I'm here!"<<endl;
	TLorentzVector ll=l1+l2,j0,j1,j2;

	}
	jetVeto->clear();
	delete jetVeto;
}


//------------------- Produce Plots DATA ---------------------------------
void QGPlotData(const char *fileIn,const char *fileOut, const char *dirOut){
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
void QGPlotMC(const char *fileIn,const char *fileOut, const char *dirOut){
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

Read A(configFile.c_str());

sscanf(A.ReadParameter("JETPT"),"%f",&JetPtCut );
sscanf(A.ReadParameter("JETDR"),"%f",&JetDRCut);
//sscanf(A.ReadParameter("LLM"),"%f",&llMCut ); //Met Cut is minimal
sscanf(A.ReadParameter("CHID2"),"%d",&CHID2 );

printf("********CUT********\n");
printf("* JetPt=%4.1f      *\n",JetPtCut);
printf("* JetDR=%4.2f      *\n",JetDRCut);
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
	if(forkNum==1)QGPlotMC( (DirMC+DY).c_str(),Form("QG_DY_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==2)QGPlotMC( (DirMC+TT).c_str(),Form("QG_TT_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==3)QGPlotMC( (DirMC+WJ).c_str(),Form("QG_WJ_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==4)QGPlotMC( (DirMC+WW).c_str(),Form("QG_WW_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==5)QGPlotMC( (DirMC+WZ).c_str(),Form("QG_WZ_%d.root",-CHID2), DirOut.c_str());
	if(forkNum==6)QGPlotMC( (DirMC+ZZ).c_str(),Form("QG_ZZ_%d.root",-CHID2), DirOut.c_str());
	
	if(CHID2==-4){if(forkNum==7)QGPlotData( (DirData+DoubleMu).c_str(),"QG_DoubleMu_4.root", DirOut.c_str());}
	if(CHID2==-1){if(forkNum==7)QGPlotData( (DirData+DoubleE).c_str(),"QG_DoubleE_1.root", DirOut.c_str());}
	if(CHID2==-2){if(forkNum==7)QGPlotData( (DirData+MuEG).c_str(),"QG_MuEG_2.root", DirOut.c_str());}


} else { //!FORK
	QGPlotMC( (DirMC+DY).c_str(),Form("QG_DY_%d.root",-CHID2), DirOut.c_str());
	QGPlotMC( (DirMC+TT).c_str(),Form("QG_TT_%d.root",-CHID2), DirOut.c_str());
	QGPlotMC( (DirMC+WJ).c_str(),Form("QG_WJ_%d.root",-CHID2), DirOut.c_str());
	QGPlotMC( (DirMC+WW).c_str(),Form("QG_WW_%d.root",-CHID2), DirOut.c_str());
	QGPlotMC( (DirMC+WZ).c_str(),Form("QG_WZ_%d.root",-CHID2), DirOut.c_str());
	QGPlotMC( (DirMC+ZZ).c_str(),Form("QG_ZZ_%d.root",-CHID2), DirOut.c_str());
	
	if(CHID2==-4){QGPlotData( (DirData+DoubleMu).c_str(),"QG_DoubleMu_4.root", DirOut.c_str());}
	if(CHID2==-1){QGPlotData( (DirData+DoubleE).c_str(),"QG_DoubleE_1.root", DirOut.c_str());}
	if(CHID2==-2){QGPlotData( (DirData+MuEG).c_str(),"QG_MuEG_2.root", DirOut.c_str());}

}
}
#endif
