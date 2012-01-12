#include "TTree.h"
#include "TFile.h"
#include <stdio.h>
#include <stdlib.h>

int removeFromDataSet(const char *fileName,const char *outFileName)
{
TFile *f=TFile::Open(fileName);if(f==NULL)return 1;

TTree *t=(TTree*)f->Get("jetTree"); if(t==NULL)return 2;

TFile *F=TFile::Open(outFileName,"RECREATE");if(F==NULL)return 3;

//t->SetBranchStatus("*",0);//smaller Tree
//t->SetBranchStatus("nCharged",1);
//t->SetBranchStatus("nNeutral",1);
//t->SetBranchStatus("eventWeight",1);
//t->SetBranchStatus("ptD",1);
//t->SetBranchStatus("rmsCand",1);
//t->SetBranchStatus("Likelihood",1);
//t->SetBranchStatus("ptJetReco",1);
//t->SetBranchStatus("etaJetReco",1);
//t->SetBranchStatus("rhoPF",1);

//PhotonJet_2ndLevelTreeW_QCD_Pt-1000to1400_pfakt5.root
//getting max pthat from dataset name
float PtHatMax,PtHatMin;
sscanf(fileName,"PhotonJet_2ndLevelTreeW_QCD_Pt-%fto%f_pfakt5.root",&PtHatMin,&PtHatMax);
char cut[1023];
sprintf(cut,"ptJetReco<%f",PtHatMax);
printf("Sel Ratio=%lf\n",t->GetEntries(cut)/float(t->GetEntries()));
TTree *T=t->CopyTree(cut);
//T->Print();
F->Write();
delete f;
delete F;
return 0;
}
