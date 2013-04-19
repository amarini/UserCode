#include <stdio.h>
#include <string>
#include "BaseAnalysis.h"
#define MAX_STR_LENGTH 1023
//==============================
char* BaseAnalysis::Read(const char*fileName,const char * parName,char*R)
{
FILE *fr=fopen(fileName,"r");
if(fr==NULL) return NULL;
if(R==NULL)R=new char[MAX_STR_LENGTH]; 
char P[MAX_STR_LENGTH];
char S[MAX_STR_LENGTH];

//leggi una linea su S
int STATUS=1;
while(STATUS!=EOF)
{
	char c='\0';
	int i=0;
	while ( (STATUS=fscanf(fr,"%c",&c))!=EOF ){
		if(c=='\0') break;
		if(c=='\n') break;
		S[i]=c;i++;
		}
	S[i]='\0';
	i=0;
	while(S[i]!='=' && S[i]!='\0'){P[i]=S[i];i++;}
	P[i]='\0';int j=0;
	if(S[i]=='=')i++;
	while(S[i]!='\0'){R[j]=S[i];i++;j++;}
	R[j]='\0';
	bool isPar=true;
	for(i=0;;i++)
		{
		if(P[i]=='\0' && parName[i]=='\0')break;
		if(P[i]=='\0' || parName[i]=='\0'){isPar=false; break;} // but not together
		if(parName[i]!=P[i]){isPar=false;break;}
		}
	
	if(isPar){fclose(fr);return R;}
}
fclose(fr);
return NULL;
}

char* BaseAnalysis::ReadMult(const char * fileName1,const char*fileName2,const char*parName){
	char* par1=Read(fileName1,parName);
	if(fileName2[0]=='\0') return par1;

	char* par2=Read(fileName2,parName);
	if(par2!=NULL)return par2;
	if(par1!=NULL)return par1;
	return NULL;
}
