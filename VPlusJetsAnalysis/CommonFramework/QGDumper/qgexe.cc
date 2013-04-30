#include "QGDumper.h"

#include <stdio.h>
int main(int argc, char**argv){

//config2 type j nJ
if(argc<3) {printf("missing parameters\nUsage: config config type j nJ\n") ; return 0;}

int type; sscanf(argv[3],"%d",&type);

int j=0; 
int nJ=1;
if(argc>5){
   sscanf(argv[4],"%d",&j);
   sscanf(argv[5],"%d",&nJ);
}

QGDumper *A=new QGDumper();
	
	printf("SetConfig\n");
	A->SetConfigFile("data/config.ini",argv[1]);
	A->SetCuts();
	A->treeName="tree_passedEvents";
	A->outFileName=Form("%s/QGDump_type%d_chid%d_%d_%d.root",A->ReadMult(argv[1],argv[2],"OUTDIR"),type,A->GetChId(),j,nJ);
	printf("BookFile\n");
	switch(type){
		case 0:	if(A->GetChId()==4)A->BookFile("DoubleMu"); //create a switch between these
			if(A->GetChId()==1)A->BookFile("DoubleE");
			if(A->GetChId()==2)A->BookFile("EMu"); break;
		case 1: A->BookFile("DY"); break;
		case 2: A->BookFile("TT"); break;
		case 3: A->BookFile("WJ"); break;
		case 4: A->BookFile("WZ"); break;
		case 5: A->BookFile("ZZ"); break;
		case 6: A->BookFile("WW"); break;
	}
	printf("OpenFiles\n");
	A->OpenFiles(); //use the booked files: will loop on the TChain

	printf("Create Tree\n");
	A->CreateHistos();

	printf("Loop\n");
	A->Loop(type,j,nJ); //type=0 data 1 MC; j , nJobs
	

	printf("Write\n");
	//A->Write(Form("%s/output_type%d_chid%d_%d_%d.root",A->ReadMult("data/config.ini",argv[1],"OUTDIR"),type,A->GetChId(),j,nJ));
	A->Write();

}
