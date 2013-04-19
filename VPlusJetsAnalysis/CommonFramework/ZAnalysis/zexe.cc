#include "ZAnalysis.h"

int main(){

ZAnalysis *A=new ZAnalysis();

	A->SetConfigFile("data/config.ini","data/config_MM.ini");
		A->BookFile("DoubleMu"); //create a switch between these
		//A->BookFile("DoubleE");
		//A->BookFile("EMu");
		//A->BookFile("DY");
		//A->BookFile("TT");
		//A->BookFile("WJ");
		//A->BookFile("WZ");
		//A->BookFile("ZZ");
		//A->BookFile("WW");
	A->OpenFiles(); //use the booked files: will loop on the TChain
	A->CreateHistos();
	A->Loop(0,1,100); //type=0 data 1 MC; j , nJobs
	
	A->SetSmearType(1,"_JES_UP");
	A->CreateHistos();
	A->Loop(0,1,100); //type=0 data 1 MC; j , nJobs

	A->SetSmearType(2,"_JES_DN");
	A->CreateHistos();
	A->Loop(0,1,100); //type=0 data 1 MC; j , nJobs

	A->Write("output_DoubleMu_1.root");

}
