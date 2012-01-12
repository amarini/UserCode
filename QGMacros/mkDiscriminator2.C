#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
/* This version wants to improve efficiency
 *
 */
using namespace std;
int getBins(double  *Bins,int nBins,double MinBin=15.0,double MaxBin=1000.,bool log=false);
int getBin(int nBins,double  *Bins,double value,double*x0=0,double*x1=0);

int mkDiscriminator2(const char *fileName="QCD_all.root",const char *variables="ptD nCharged nNeutral rmsCand",const char *treeName="jetTree")
{
double PtBins[1023];
double RhoBins[1023];
int nRhoBins=20;
int nPtBins=20;
getBins(PtBins,nPtBins,15,1000.,true);
getBins(RhoBins,nRhoBins,0,20.,false);
map< string, TH1F *> plots;
TFile *f=TFile::Open(fileName);
TTree *t=(TTree*)f->Get(treeName);
TFile *F=TFile::Open("outFile.root","RECREATE");

//general things of interest
float ptJetReco;
float etaJetReco;
int   pdgIdPart;
void  *Variable;
float eventWeight;
float rhoPF;

t->SetBranchAddress("ptJetReco",&ptJetReco);
t->SetBranchAddress("etaJetReco",&etaJetReco);
t->SetBranchAddress("pdgIdPart",&pdgIdPart);
t->SetBranchAddress("eventWeight",&eventWeight);
t->SetBranchAddress("rhoPF",&rhoPF);

char str[1023];
char cut[1023];
char plotName[1023];
char VarName[1023];
const char *VariablesPointer=variables; int n;
fprintf(stderr,"Beginning var loops\n");
while(sscanf(VariablesPointer,"%s%n",VarName,&n)==1)
	{
	fprintf(stderr,"Variable=%s\n",VarName);
	VariablesPointer+=n;
	int nBinsX=100;
	float xMin=0,xMax=1;
	switch(VarName[1])
	{
	case 't': nBinsX=50;xMin=0;xMax=1.0;Variable=new float;t->SetBranchAddress(VarName,Variable);break;//ptD
	case 'C': nBinsX=101;xMin=-.5;xMax=100.5;Variable=new int  ;t->SetBranchAddress(VarName,Variable);break;//nCharged
	case 'N': nBinsX=101;xMin=-.5;xMax=100.5;Variable=new int  ;t->SetBranchAddress(VarName,Variable);break;//nNeutral
	case 'm': nBinsX=100;xMin=0;xMax=1.0;Variable=new int  ;t->SetBranchAddress(VarName,Variable);break;//rmsCand
	default: break;
	}
	//for each bins in pt
	for(int p=0;p<nPtBins;p++)
	{
	//create the directory
	F->cd();
	sprintf(str,"rhoBins_pt%.0lf_%.0lf",ceil(PtBins[p]),ceil(PtBins[p+1]));
	F->mkdir(str);
	F->cd(str);
	for(int r=0;r<nRhoBins;r++)
	{
	sprintf(plotName,"%s_gluon_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));//construction of plot name
	plots[string(plotName)]=new TH1F(plotName,plotName,nBinsX,xMin,xMax);
	sprintf(str,"%s>>%s",VarName,plotName);
	sprintf(cut,"eventWeight*(pdgIdPart==21 && %lf< ptJetReco && ptJetReco<%lf && %lf< rhoPF && rhoPF<%lf && abs(etaJetReco)<2.4)",PtBins[p],PtBins[p+1],RhoBins[r],RhoBins[r+1]);
	//quark
	sprintf(plotName,"%s_quark_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));//construction of plot name

	plots[string(plotName)]=new TH1F(plotName,plotName,nBinsX,xMin,xMax);
	sprintf(str,"%s>>%s",VarName,plotName);
	//quark = uds not cb
	sprintf(cut,"eventWeight*(abs(pdgIdPart)<4 && %lf< ptJetReco && ptJetReco<%lf && %lf< rhoPF && rhoPF<%lf && abs(etaJetReco)<2.4 )",PtBins[p],PtBins[p+1],RhoBins[r],RhoBins[r+1]);
	}//loop on rho Bins
	}//loop on pt bins
	//tree loop
	for(long long int entry=0;entry<t->GetEntries();entry++)
	{
		t->GetEntry(entry);
		double ptBin0,ptBin1,rhoBin0,rhoBin1;
		if(getBin(nPtBins,PtBins,ptJetReco,&ptBin0,&ptBin1)<0)continue;
		if(getBin(nRhoBins,RhoBins,rhoPF,&rhoBin0,&rhoBin1)<0)continue;
		//construct plotname
		if(pdgIdPart==21)sprintf(plotName,"%s_gluon_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(ptBin0),ceil(ptBin1),floor(rhoBin0));
		else if((-4<pdgIdPart)&&(pdgIdPart)<4)sprintf(plotName,"%s_quark_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(ptBin0),ceil(ptBin1),floor(rhoBin0));
		else continue;
		//selection
		if(etaJetReco>2.0)continue;
		if(etaJetReco<-2.0)continue;
		//Fill
		//plots[PlotName]->Fill(*(int)Variable,eventWeight);
		switch(VarName[1])
		{
		case 't':plots[string(plotName)]->Fill(*(float*)Variable,eventWeight); break;//ptD
		case 'C':plots[string(plotName)]->Fill(*(int*)Variable,eventWeight); break;//nCharged
		case 'N':plots[string(plotName)]->Fill(*(int*)Variable,eventWeight); break;//nNeutral
		case 'm':plots[string(plotName)]->Fill(*(float*)Variable,eventWeight); break;//rmsCand
		default: break;
		}
	}//loop on entries	

	}//loop on variables names
//Write plots
map< string, TH1F *>::iterator plots_iterator;
for(plots_iterator=plots.begin();plots_iterator!=plots.end();plots_iterator++)
	{
	float Pt0,Pt1,Rho0;
	char Dir[1023],var[1023],tmp[1023],pdg[1023];
	fprintf(stderr,"Writing file %s\n",plots_iterator->first.c_str());
	sscanf(plots_iterator->first.c_str(),"%s",tmp);
	//sscanf works properly with %s with space on newline ...
	for(int i=0;;i++)if(tmp[i]=='_')tmp[i]=' '; else if(tmp[i]=='\0') break;
	sscanf(tmp,"%s %s pt%f %f rho%f",pdg,var,&Pt0,&Pt1,&Rho0);//no ceil	
	sprintf(Dir,"rhoBins_pt%.0f_%.0f",Pt0,Pt1);
	fprintf(stderr,"  Dir: %s\n",Dir);
	F->cd(Dir);
	plots_iterator->second->Write();
	}
return 0;
}

int getBins(double  *Bins,int nBins,double MinBin,double MaxBin,bool log)
{
double incr;
if(log)
	{
	incr=TMath::Power(MaxBin/MinBin,1.0/double(nBins));
	Bins[0]=MinBin;
	Bins[nBins]=MaxBin;
	for(int i=1;i<nBins;i++)
		Bins[i]=Bins[i-1]*incr;
	}
else
	{
	incr=(MaxBin-MinBin)/nBins;
	Bins[0]=MinBin;
	Bins[nBins+1]=MaxBin;
	for(int i=1; i<nBins+1;i++)
		Bins[i]=Bins[i-1]+incr;
	}
return 0;
}
int getBin(int nBins,double  Bins[],double value,double *x0,double *x1)
{
int R=0;
//int nBins=sizeof(Bins)/sizeof(double);//?
if(value <Bins[0])return -1;
if(value >Bins[nBins])return -1;
for(R=0;R<nBins;R++)
	{
	if(Bins[R]>value)break;	
	}
R--;
if(x0) *x0=Bins[R];
if(x1) *x1=Bins[R+1];
return R;	
}
