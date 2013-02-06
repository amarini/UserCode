#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1I.h"
//unit definition of femtobarn vs picobarn: it may be useful to have plot normalized to 1fb
#define FB 1000.

#include <stdio.h>
#include <stdlib.h>

int AddEventWeight(	const char * FileName="ntuple.root", //File name
			const char* Directory="zjetgen", //Directory name in the root file 
			const double CrossSection=1300, // Cross Section of the Sample
			const char* TreeName="events", // Tree Name in the directory Chosen
			const char* HistoName="WEvents",  // histo name in the directory Chosen
			const float kFactor=1.0
		   	)
{
	//declare a temporary string variable
	char str[255];	
	//open file
	TFile *f=TFile::Open(FileName,"UPDATE");
	if(f==NULL){
		fprintf(stderr,"%s: No such file or directory\n",FileName);
		return 1;
		}
	//cd to directory
	f->cd(Directory);
	//Minimum check on Cross Section
	if(CrossSection<0.0){
		fprintf(stderr,"Cross Section (%lf) <0\n",CrossSection);
		return 2;
		}
	//Constructing TreeName and open it
	sprintf(str,"%s/%s",Directory,TreeName);
	TTree *t=(TTree*)f->Get(str);
	if(t==NULL){
		fprintf(stderr,"%s: No such tree\n",TreeName);
		return 3;
		}
	//Constructing HistoName and open it 
	sprintf(str,"%s/%s",Directory,HistoName);
	TH1I *h=(TH1I*)f->Get(str);
	if(h==NULL){
		fprintf(stderr,"%s: No such histo\n",HistoName);
		return 4;
		}
	//Get the Number of event Processed 
	double NumberEventProcessed=h->Integral();
	//Compute the event Weight as function of cross section and number of event processed 
	double eventWeight=CrossSection*FB/NumberEventProcessed;
	//Creating an empty branch in the tree
	TBranch *b=t->Branch("eventWeight",&eventWeight,"eventWeight/D");
	//setting the branch for the mcWeight
	float mcWeight=1.0;
	t->SetBranchAddress("mcWeight",&mcWeight);
	//Getting the Number of entries in the tree
	long long int NumberEntries=t->GetEntries();
	//looping on the entries in order to add the correct number of entries in the branch
	for(long long int i=0;i<NumberEntries;i++){
		t->GetEntry(i);
		eventWeight=double((CrossSection*FB/NumberEventProcessed)) * double(mcWeight) *kFactor;
		b->Fill();
		}
	//Write the Tree (With OverWrite Option)
	t->Write("",TObject::kOverwrite);
	//Close the file
	f->Close();
	//Print a message on stdout
	printf("Done\n");
	return 0;
}

