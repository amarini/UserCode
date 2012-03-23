#include "PtMeans.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "PtBins.C"

#define MAXBINS 100
Means::Means(){
PtMeans =new double[MAXBINS];
RhoMeans=new double[MAXBINS];
}

Means::~Means(){
delete[] PtMeans;
delete[] RhoMeans;
}

int Means::ReadTxt(const char *fileName)
{
FILE *fr=fopen(fileName,"r");
if(fr==NULL)return 1;
char buffer[1023];
fscanf(fr,"%[^\n]\n", buffer);//TO CHECK
int i=0;while(sscanf(buffer,"%lf",&PtMeans[i]))i++;
fscanf(fr,"%[^\n]\n", buffer);
i=0;while(sscanf(buffer,"%lf",&RhoMeans[i]))i++;
return 0;
}

int Means::WriteTxt(const char *fileName,const char*outFileName,const char *treeName)
{
//this file is supposted to be omog
TFile *f=TFile::Open(fileName);
TTree *t=(TTree*)f->Get(treeName);
FILE *fw=fopen(outFileName,"w");
if(fw==NULL) return 1;
double PtBins[MAXBINS];
double RhoBins[MAXBINS];
//actual bin composition -- TODO Write this in BINS
getBins_int(18,PtBins,20,1000,true);
PtBins[18]=3500;
getBins_int(21,RhoBins,0,20,false);
//Compute PtMeans Rho Integrated
for( int i=0; i<18;i++){
		char selection[1023];
		sprintf(selection,"eventWeight*(%f<ptJet0 && ptJet0<%f)",PtBins[i],PtBins[i+1]);
		TH1F *tmp=new TH1F("tmp","tmp",3500,0,3500);
		t->Draw("ptJet0>>tmp",selection);
		PtMeans[i]=tmp->GetMean();
		PtMeans[i+1]=-1;
		delete tmp;
		}
//Compute RhoMeans Pt Integrated
for(int  i=0; i<20;i++){
		char selection[1023];
		sprintf(selection,"eventWeight*(%f<rhoPF && rhoPF<%f)",RhoBins[i],RhoBins[i+1]);
		TH1F *tmp=new TH1F("tmp","tmp",300,0,30);
		t->Draw("rhoPF>>tmp",selection);
		RhoMeans[i]=tmp->GetMean();
		RhoMeans[i+1]=-1;
		delete tmp;
		}
//Write on file
int i;
i=0;while(PtMeans[i]>=0.0){fprintf(fw,"%.5f ",PtMeans[i]);i++;}
nPtMeans=i+1;
fprintf(fw,"\n");
i=0;while(RhoMeans[i]>=0.0){fprintf(fw,"%.5f ",RhoMeans[i]);i++;}
nRhoMeans=i+1;
}
