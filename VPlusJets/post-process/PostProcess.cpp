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
			}	
		}
	if(m==1) return R;
	else if( m==0) return CrossSection::noMatch;
	else if( m>1) return CrossSection::multipleMatch;
	}

#include "AddEventWeight.C"
int Usage(const char *progName){
	printf("Usage:\n");
	printf("      %s fileName\n",progName);
	return 0;
}

int main(int argc, char**argv){
	if(argc<2){return Usage(argv[0]);}
	if(argc==2){
		CrossSection A;
		A.ReadTxtFile("xSec.ini");
		double xSec=A.xSection(argv[1]);
		if (xSec==CrossSection::noMatch){printf("No match Cross Section in the Database\n");return -1;}
		if (xSec==CrossSection::multipleMatch){printf("Multiple matches in Cross Section in the Database\n");return -1;}

		return AddEventWeight(argv[1],"accepted",xSec,"events","WEvents",1.0);
		}

}
