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
			const char * What="L C N P R",//what I will add F:LikelihoodFit
			const char* TreeName="jetTree", // Tree Name in the directory Chosen
			const char *QGFile=""
		   	)
{
	//reading what
	const char *pointer=What;
	int n;char a;
	bool L=false,C=false,N=false,P=false,R=false,F=false;
	while(sscanf(pointer,"%c %n",&a,&n)==1)
		{
		pointer+=n;
		switch (a){
		case 'L': L=true;fprintf(stderr,"L ");break;
		case 'C': C=true;fprintf(stderr,"C ");break;
		case 'N': N=true;fprintf(stderr,"N ");break;
		case 'P': P=true;fprintf(stderr,"P ");break;
		case 'R': R=true;fprintf(stderr,"R ");break;
		case 'F': F=true;fprintf(stderr,"F ");break;
		}
		};
	fprintf(stderr,"\n");	
	//
	fprintf(stderr,"creating the qgDiscr\n");
	QGLikelihoodCalculator *qglikeli;
	qglikeli = new QGLikelihoodCalculator("/Users/andreamarini/Documents/QGDiscriminator/FitResults3/nCharged_2ndLevel.txt","/Users/andreamarini/Documents/QGDiscriminator/FitResults3/nNeutral_2ndLevel.txt","/Users/andreamarini/Documents/QGDiscriminator/FitResults3/ptD.txt");
	
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
	float Likelihood;
	float LikelihoodFit;
	float ptD;
	float rmsCand;
	int nCharged;
	int nNeutral;
	//Creating an empty branch in the tree
	fprintf(stderr,"Creating the branches and setting the address\n");
	TBranch *b0;if(L)b0=t->Branch("Likelihood",&Likelihood,"Likelihood/F");
	TBranch *b1;if(P)b1=t->Branch("ptD",&ptD,"ptD/F"); //used ptD
	TBranch *b2;if(C)b2=t->Branch("nCharged",&nCharged,"nCharged/I"); //used nCharged
	TBranch *b3;if(N)b3=t->Branch("nNeutral",&nNeutral,"nNeutral/I"); //used nNeutral
	TBranch *b4;if(R)b4=t->Branch("rmsCand",&rmsCand,"rmsCand/F"); //used rmsCand
	TBranch *b5;if(F)b5=t->Branch("QGFit2",&LikelihoodFit,"QGFit2/F"); //used LikelihoodFit
	
	//Getting the Number of entries in the tree
	long long int NumberEntries=t->GetEntries();
	//Setting Branch Address
	float ptJetReco            ;
	float etaJetReco           ;
	float rhoPF                ;
	float ptDJetReco           ;
	float rmsCandJetReco           ;
	int   nNeutralHadronsReco  ;
	int   nPhotonsReco         ;
	int   nMuonsReco           ;
	int   nElectronsReco       ;
	int   nTracksReco          ;
	t->SetBranchAddress("ptJet0",&ptJetReco);
	t->SetBranchAddress("rmsCandJet0",&rmsCandJetReco);
	t->SetBranchAddress("etaJet0",&etaJetReco);
	t->SetBranchAddress("rhoPF",&rhoPF);
	t->SetBranchAddress("ptDJet0",&ptDJetReco);
	if(N)t->SetBranchAddress("nNeutralHadronsReco",&nNeutralHadronsReco);
	if(N)t->SetBranchAddress("nPhotonsReco",&nPhotonsReco);
	if(C)t->SetBranchAddress("nMuonsReco",&nMuonsReco);
	if(C)t->SetBranchAddress("nElectronsReco",&nElectronsReco);
	if(C)t->SetBranchAddress("nTracksReco",&nTracksReco);
	if(!C)t->SetBranchAddress("nChargedJet0",&nCharged);
	if(!N)t->SetBranchAddress("nNeutralJet0",&nNeutral);
	//looping on the entries in order to add the correct number of entries in the branch
	fprintf(stderr,"Beginning the loop over the entries\n");
	for(long long int i=0;i<NumberEntries;i++){
		t->GetEntry(i);
		if(P)ptD=ptDJetReco;
		if(R)rmsCand=rmsCandJetReco;
		if(C)nCharged=nMuonsReco+nElectronsReco+nTracksReco;
		if(N)nNeutral=nNeutralHadronsReco+nPhotonsReco;
		if(L)Likelihood   =qglikeli->computeQGLikelihoodPU(ptJetReco,rhoPF,nCharged,nNeutral,ptDJetReco);
		if(F)LikelihoodFit=qglikeli->computeQGLikelihoodPU(ptJetReco,rhoPF,nCharged,nNeutral,ptDJetReco);

		if((i&131071)==1)fprintf(stderr,"entry %lld of %lld: Likelihood=%f,ptJet=%f,rho=%f, LikelihoodFit=%f\n",i,NumberEntries,Likelihood,ptJetReco,rhoPF,LikelihoodFit);
		//fprintf(stderr,"entry %lld of %lld: Likelihood=%f,ptJet=%f,rho=%f, LikelihoodFit=%f\n",i,NumberEntries,Likelihood,ptJetReco,rhoPF,LikelihoodFit);
		
		if(L)b0->Fill();
		if(P)b1->Fill();
		if(C)b2->Fill();
		if(N)b3->Fill();
		if(R)b4->Fill();
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
