#include <cstdio>
#include <stdlib.h>
#include <string>
#include <map>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1I.h"

using namespace std;

class CrossSection{
public:
	int ReadTxtFile(const char*fileName);
	const static int noFile=-1;
	const static double noMatch=-2;
	const static double multipleMatch=-3;
	double xSection(string match);
private:
	map<string,double> xSec;	
};
int CrossSection::ReadTxtFile(const char*fileName)
	{
	FILE *fr=fopen(fileName,"r");
	if(fr==NULL) return CrossSection::noFile;
	char key[1023];
	double number;
	while(fscanf(fr,"%s %lf",key,&number)!=EOF){
		xSec[string(key)]=number;
		}
	return 0;
	}
double CrossSection::xSection(string match){
	float R=0;
	int m=0;
	for(map<string,double>::iterator it=xSec.begin();it!=xSec.end();it++)
		{
		if(match.find(it->first) != string::npos){ // I want to match what I put in the database and not the fileName
			R=it->second;
			m++;
			//DEBUG
			fprintf(stderr,"MATCH=%s\n",it->first.c_str());
			}	
		}
	if(m==1) return R;
	else if( m==0) return CrossSection::noMatch;
	else if( m>1) return CrossSection::multipleMatch;
	}

#include "AddEventWeight.C"
#include "AddPUWeight.C"

int Usage(const char *progName){
	printf("Usage:\n");
	printf("      %s fileName [pu]\n",progName);
	printf("             - pu 0: add Cross Section Weights  And PUWeight. Pileupfile is pileup/all.json.pileup.root\n");
	printf("             - pu 1: add PUWeight+1s. Pileupfile is pileup/all.json.pileup_UP.root\n");
	printf("             - pu -1: add PUWeight-1s. Pileupfile is pileup/all.json.pileup_DN.root\n");
	printf("             - pu 10: add PUWeight \n");
	return 0;
}

int main(int argc, char**argv){
	if(argc<2){return Usage(argv[0]);}
	if(argc==2){
		printf("2\n");//DEBUG
		CrossSection A;
			size_t found;
			found=string(argv[0]).find_last_of("/\\");	
			string folder=string(argv[0]).substr(0,found);
		A.ReadTxtFile( (folder+ "/xSec.ini").c_str());
		double xSec=A.xSection(argv[1]);
		if (xSec==CrossSection::noMatch){printf("No match Cross Section in the Database\n");return -1;}
		if (xSec==CrossSection::multipleMatch){printf("Multiple matches in Cross Section in the Database\n");return -1;}
		printf("xSec=%.1lf\n",xSec);
		return AddEventWeight(argv[1],"accepted",xSec,"events","WEvents",1.0);
		}
	if(argc==3){
		int pu;
		printf("3\n");//DEBUG
		sscanf(argv[2],"%d",&pu);
		if(pu==0){
			CrossSection A;
				size_t found;
				found=string(argv[0]).find_last_of("/\\");	
				string folder=string(argv[0]).substr(0,found);
			A.ReadTxtFile( (folder+ "/xSec.ini").c_str());
			double xSec=A.xSection(argv[1]);
			if (xSec==CrossSection::noMatch){printf("No match Cross Section in the Database\n");return -1;}
			if (xSec==CrossSection::multipleMatch){printf("Multiple matches in Cross Section in the Database\n");return -1;}
			printf("xSec=%.1lf\n",xSec);
			AddEventWeight(argv[1],"accepted",xSec,"events","WEvents",1.0);
			AddPUWeight(argv[1],"mcPU");
		}else if(pu==10){
			AddPUWeight(argv[1],"mcPU");
		}else if(pu ==1){
			AddPUWeight(argv[1],"mcPU","accepted","events",1);
		}else if(pu==-1){
			AddPUWeight(argv[1],"mcPU","accepted","events",-1);
			}
			
		
		}

}
