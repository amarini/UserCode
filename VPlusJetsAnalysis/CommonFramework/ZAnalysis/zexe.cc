#include "ZAnalysis.h"

#include <stdio.h>
int main(){

ZAnalysis *A=new ZAnalysis();
	
	printf("SetConfig\n");
	A->SetConfigFile("data/config.ini","data/config_MM.ini");
	A->SetCuts();
	printf("BookFile\n");
		A->BookFile("DoubleMu"); //create a switch between these
		//A->BookFile("DoubleE");
		//A->BookFile("EMu");
		//A->BookFile("DY");
		//A->BookFile("TT");
		//A->BookFile("WJ");
		//A->BookFile("WZ");
		//A->BookFile("ZZ");
		//A->BookFile("WW");
	printf("OpenFiles\n");
	A->OpenFiles(); //use the booked files: will loop on the TChain
	printf("Create Histos\n");
	A->CreateHistos();
	printf("Loop\n");
	A->Loop(0,1,100); //type=0 data 1 MC; j , nJobs
	
	printf("Smear JES UP\n");
	A->SetSmearType(1,"_JES_UP");
	printf("--> Histos\n");
	A->CreateHistos();
	printf("--> Loop\n");
	A->Loop(0,1,100); //type=0 data 1 MC; j , nJobs

	printf("Smear JES DN\n");
	A->SetSmearType(2,"_JES_DN");
	printf("--> Histos\n");
	A->CreateHistos();
	printf("--> Loop\n");
	A->Loop(0,1,100); //type=0 data 1 MC; j , nJobs

	printf("Write\n");
	A->Write("output_DoubleMu_1.root");

}
