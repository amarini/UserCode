#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
//unit definition of femtobarn vs picobarn: it may be useful to have plot normalized to 1fb

#include <stdio.h>
#include <stdlib.h>
#include <vector>

//#include "QGLikelihood/QGLikelihoodCalculator.C"

#include "QGLikelihoodCalculator.C"

using namespace std;

int AddLikelihoodFit(	const char * FileName="ntuple.root", //File name
			const char * What="F",//what I will add F:LikelihoodFit
			const char* TreeName="jetTree", // Tree Name in the directory Chosen
			const char *QGFile=""
		   	)
{
	//reading what
	const char *pointer=What;
	int n;char a;
	bool F=false;
	while(sscanf(pointer,"%c %n",&a,&n)==1)
		{
		pointer+=n;
		switch (a){
		case 'F': F=true;fprintf(stderr,"F ");break;
		}
		};
	fprintf(stderr,"\n");	
	//
	fprintf(stderr,"creating the qgDiscr\n");
	QGLikelihoodCalculator *qglikeli;
	if(QGFile[0]=='\0')	qglikeli = new QGLikelihoodCalculator();
	else 			qglikeli = new QGLikelihoodCalculator(QGFile);
	
	//declare a temporary string variable
	char str[255];	
	//open file
	fprintf(stderr,"Opening the file and gettin the tree\n");
	TFile *f=TFile::Open(FileName,"UPDATE");
	if(f==NULL){
		fprintf(stderr,"%s: No such file or directory\n",FileName);
		return 1;
		}
	//Constructing TreeName and open it
	TTree *t=(TTree*)f->Get(TreeName);
	if(t==NULL){
		fprintf(stderr,"%s: No such tree\n",TreeName);
		return 2;
		}
	float LikelihoodFit;
	float ptD;
	int nCharged;
	int nNeutral;
	//Creating an empty branch in the tree
	fprintf(stderr,"Creating the branches and setting the address\n");
	TBranch *b5;if(F)b5=t->Branch("QGFit2",&LikelihoodFit,"QGFit2/F"); //used LikelihoodFit
	
	//Getting the Number of entries in the tree
	long long int NumberEntries=t->GetEntries();
	//Setting Branch Address
	float rhoPF                ;
	float ptJetReco                ;
	t->SetBranchAddress("rhoPF",&rhoPF);
	t->SetBranchAddress("ptJetReco",&ptJetReco);
	t->SetBranchAddress("ptD",&ptD);
	t->SetBranchAddress("nCharged",&nCharged);
	t->SetBranchAddress("nNeutral",&nNeutral);
	//t->SetBranchAddress("ptD",&ptD);
	//looping on the entries in order to add the correct number of entries in the branch
	fprintf(stderr,"Beginning the loop over the entries\n");
	for(long long int i=0;i<NumberEntries;i++){
		t->GetEntry(i);
		if(F)LikelihoodFit=qglikeli->computeQGLikelihoodPU(ptJetReco,rhoPF,nCharged,nNeutral,ptD);

		if((i&131071)==1)fprintf(stderr,"entry %lld of %lld: LikelihoodFit=%f,ptJet=%f,rho=%f\n",i,NumberEntries,LikelihoodFit,ptJetReco,rhoPF);
		
		if(F)b5->Fill();
		}
	//Write the Tree (With OverWrite Option)
	t->Write("",TObject::kOverwrite);
	//Close the file
	f->Close();
	//Print a message on stdout
	printf("Done\n");
	return 0;
}
