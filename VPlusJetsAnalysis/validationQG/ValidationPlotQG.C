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

float JetPtCut=30;
float JetDRCut=0.4;
float llMCut=20;
int CHID2=-4;
#define MAXBINS 30
float PtBins[MAXBINS];
float RhoBins[MAXBINS];


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

//rho
histos["rho"]		= new TH1F("rho","rho;#rho;events",50,0,50);
histos["rho25"]		= new TH1F("rho25","rho25;#rho(2.5); events",50,0,50);

//jets
histos["jet0Pt"]	= new TH1F("jet0Pt","jet0Pt;P_{T}^{1st jet};events",50,50,150);
for(int pt_bin=0;pt_bin+1 <MAXBINS && PtBins[pt_bin+1]>=0;pt_bin++)
for(int rho_bin=0;rho_bin+1 <MAXBINS && RhoBins[rho_bin+1]>=0;rho_bin++)
	{
	string name=Form("QGL_jet0Pt_bin%d_rho_bin%d",pt_bin,rho_bin);
	histos[name.c_str()]	= new TH1F(name.c_str(),"QGL;QGL;events",100,-1.0,1.00001);
	name=Form("QGLMLP_jet0Pt_bin%d_rho_bin%d",pt_bin,rho_bin);
	histos[name.c_str()]	= new TH1F(name.c_str(),"QGLMLP;QGLMLP;events",100,-1.0,1.00001);
	}

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
vector<float> *jetQGLMLP=NULL	;t->SetBranchAddress("jetQGLMLP",&jetQGLMLP);
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

	//llM selection
	if(!(fabs(llM-91)<llMCut))continue;
	
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
	
	//Find Pt/Rho Bin
	int pt_bin=-1;for(int i=0;i<MAXBINS && PtBins[i+1]>=0;i++){if( PtBins[i]<(*jetPt)[jet0_BS] && (*jetPt)[jet0_BS]<PtBins[i+1]){pt_bin=i;break;}}
	int rho_bin=-1;for(int i=0;i<MAXBINS && RhoBins[i+1]>=0;i++){if( RhoBins[i]<rho && rho<RhoBins[i+1]){rho_bin=i;break;}}
	
	string name=Form("QGL_jet0Pt_bin%d_rho_bin%d",pt_bin,rho_bin);
	if(jet0_BS>=0)histos[name]->Fill( (*jetQGL)[jet0_BS],weight);
	name=Form("QGLMLP_jet0Pt_bin%d_rho_bin%d",pt_bin,rho_bin);
	if(jet0_BS>=0)histos[name]->Fill( (*jetQGLMLP)[jet0_BS],weight);

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

string PtBins_str=A.ReadParameter("PTBINS");
string RhoBins_str=A.ReadParameter("RHOBINS");
	const char *ptr=PtBins_str.c_str();int n_str,pt_bin=0;
	while(sscanf(ptr,"%f%n",&PtBins[pt_bin],&n_str)>=1){ptr+=n_str;pt_bin++;}
	ptr=RhoBins_str.c_str();int rho_bin=0;
	while(sscanf(ptr,"%f%n",&RhoBins[rho_bin],&n_str)>=1){ptr+=n_str;rho_bin++;}
	PtBins[pt_bin]=-1.;
	RhoBins[rho_bin]=-1.;

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
