#include <iostream>
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include <math.h>

#ifndef COMP_SYST_SIMPLE_H
#define COMP_SYST_SIMPLE_H


#define MAX_STR_LEN 1023

class range{
	public:
		range(float pt0,float pt1, float rho0, float rho1,const char eta0){
			ptmin=pt0;ptmax=pt1;rhomin=rho0;rhomax=rho1;eta=eta0;
		};
		range(){range(0,0,0,0,'C');}
		float ptmin;
		float ptmax;
		float rhomin;
		float rhomax;
		char eta;
};

class histo{
	public:
		histo(){N=0;lmin=-1;lmax=-1;};
		int N;
		float lmin;
		float lmax;
		std::vector<float> values;
	};

const bool operator<(const range &x,const range &y){
	if(x.ptmin !=y.ptmin) return x.ptmin<y.ptmin;
	else if(x.ptmax != y.ptmax) return x.ptmax<y.ptmax;
	else if(x.rhomin !=y.rhomin) return x.rhomin<y.rhomin;
	else if(x.rhomax != y.rhomax) return x.rhomax<y.rhomax;
	else return x.eta < y.eta;
	//else return false;
	}

class ComputeSystSimple{
public:
	//Constructor - Destructor
	ComputeSystSimple(std::string fileName="SystDatabase.txt");
	~ComputeSystSimple();
		//l is the likelihood. 
	float GetSystSimple(float pt, float rho,float eta, float l);
	float GetSystSimple(float pt, float rho,const char eta, float l);
	void SetDebug(int d){debug_=d;}
	void Print();
protected:
	char GetChar(FILE *fr);
	int debug_;
	std::map< range,float> alpha;
	std::map< range,histo> eff;
};

char ComputeSystSimple::GetChar(FILE *fr)
{
fpos_t pos;//position in the file
char c;
//get position of the line to be analized;
fgetpos(fr,&pos);
c=fgetc(fr);
fsetpos(fr,&pos); //moving back to the beginning of the line
return c;
}

ComputeSystSimple::ComputeSystSimple(std::string fileName){
//DEBUG
 SetDebug(0);
//OPEN TXT FILE 
	if(debug_>0)printf("Opening file %s\n",fileName.c_str());
 FILE *fr=fopen(fileName.c_str(),"r");	
 if(fr==NULL) {
	fprintf(stderr,"File %s does not exist\n",fileName.c_str());
	//throw nofile;
	}
 char c; 
 char str[MAX_STR_LEN]; // max length of the line
int loop=0;
 while ( (c=GetChar(fr)) != EOF)
	{
		if(debug_>2)printf("Readed Character ==%c==%d==\n",c,c);
	//ungetc(c,fr);
	fgets(str,MAX_STR_LEN,fr);
		if(debug_>2)printf("Readed String ==%s==\n",str);
	if(c=='#') continue;
	if(c=='!') {printf("%s",str);continue;}
	if(c=='\n') continue;
	if(c=='[') {loop++;continue;}
	if(debug_>2)printf("Loop ==%d==\n",loop);
	switch(loop){
		case 1:
			{ 
			float pt0,pt1,rho0,rho1,a;
			char eta;
			if(sscanf(str,"%f %f %f %f %c %f",&pt0,&pt1,&rho0,&rho1,&eta,&a) >=0  ) {
			alpha[range(pt0,pt1,rho0,rho1,eta)]=a;
			if(debug_>1)printf("Setting alpha %.0f %.0f %.0f %.0f %c %.1f\n",pt0,pt1,rho0,rho1,eta,a);
			}
			else fprintf(stderr,"Error reading line: %s\n--- Undesired lines should begin with #\n",str);
			break;
			}
		case 2:
			{
			float pt0,pt1,rho0,rho1,lmin,lmax,x;
			char eta;
			int N,n;
			char *str_ptr;
			if(sscanf(str,"%f %f %f %f %c %d %f %f%n",&pt0,&pt1,&rho0,&rho1,&eta,&N,&lmin,&lmax,&n) >=0  ) {
			str_ptr=str+n;
			eff[range(pt0,pt1,rho0,rho1,eta)].lmin=lmin;
			eff[range(pt0,pt1,rho0,rho1,eta)].lmax=lmax;
			eff[range(pt0,pt1,rho0,rho1,eta)].N=N;
			if(debug_>1)printf("Setting eff %.0f %.0f %.0f %.0f %d %.2f %.2f\n",pt0,pt1,rho0,rho1,N,lmin,lmax);
			for(int i=0;i<N;i++)
				{
				x=-1;
				sscanf(str_ptr,"%f%n",&x,&n);
				str_ptr+=n;
				eff[range(pt0,pt1,rho0,rho1,eta)].values.push_back(x);
				}
			}
			break;
			}
		default:	
		case 0: break;
		}
	if(loop>2)break;
	}
 
}
ComputeSystSimple::~ComputeSystSimple(){
}

float ComputeSystSimple::GetSystSimple(float pt, float rho,float eta, float l){
if(eta<2.5)return GetSystSimple(pt,rho,'C',l);
if(eta>=2.5)return GetSystSimple(pt,rho,'F',l);
}

float ComputeSystSimple::GetSystSimple(float pt, float rho,const char eta, float l){
	//find the range in pt
	//find the eff cut 
	float x=-999;
	bool findX=false;
	
	for(std::map< range,histo>::iterator it=eff.begin();it!=eff.end();it++){
		if(pt<it->first.ptmax)
		if(pt>it->first.ptmin)
		if(rho<it->first.rhomax)
		if(rho>it->first.rhomin)
		if(eta==it->first.eta)
			{
			int N=it->second.N;
			float lmin=it->second.lmin;
			float lmax=it->second.lmax;
			int n = floor( (l-lmin)*N/(lmax-lmin) );
			if(debug_>0)printf("Readed values %d %.1f %.1f n=%d \n",N,lmin,lmax,n);
			if (n<int(it->second.values.size()))
				{
				x=it->second.values.at(n);
				findX=true;
				}
			else {fprintf(stderr,"Error in the l value - check your config file [lmin,lmax)\n");}
			}
	}
	
	//find the range in pt
	//find alpha
	float a=-1;
	bool findA=false;
	for(std::map< range,float>::iterator it=alpha.begin();it!=alpha.end();it++){
		if(pt<it->first.ptmax)
		if(pt>it->first.ptmin)
		if(rho<it->first.rhomax)
		if(rho>it->first.rhomin)
		if(eta==it->first.eta)
			a=it->second;	
			findA=true;
	}
	//fprintf(stderr,"--x=%.2f a0=%.2f a=%.2f\n",x,a,4*a*x*(1-x));//DEBUG
	if(findX&&findA)
		return 4*a*x*(1-x);
	else return -1;
}

void ComputeSystSimple::Print(){
	printf("[alpha]\n");	
	for(std::map< range,float>::iterator it=alpha.begin();it!=alpha.end();it++){
		printf("%.0f %.0f %.0f %.0f %c %.5f\n",it->first.ptmin,it->first.ptmax,it->first.rhomin,it->first.rhomax,it->first.eta,it->second);
		}
		
	}
#endif 
//COMP_SYST_SIMPLE_H 

