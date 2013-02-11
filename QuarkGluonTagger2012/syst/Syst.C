#include "TMath.h"
#include <string>

inline float function(float x, float a ,float b)
{
//TMath::ATan( %f * TMath::Tan(TMath::Pi()*%s-TMath::Pi()/2.) + %f)/TMath::Pi() +0.5 
using namespace TMath;
//return ATan( a*Tan(Pi()*x-Pi()/2.0) + b)/Pi() +0.5;
return (ATan( a*Tan( (Pi()*x)-(Pi()/2.0)) + b)/Pi()) + 0.5;
}

float Syst(const char *tagger,float pt, float rho, float eta,float value )//"QGL","MLP"
{
using namespace std;
int PtBin=-999;
int RhoBin=-999;
int EtaBin=-999;
int tag=-999;

if(string(tagger)==string("QGL")) tag=1;
if(string(tagger)==string("MLP")) tag=2;

if( rho <=15 ) RhoBin=1;
if( rho >15  ) RhoBin=2;

if( eta <=2.5) {EtaBin=1;
	if( pt <50  ) PtBin=1;
	if( 50<= pt && pt <120 ) PtBin=2;
	if( 120<= pt ) PtBin=3;
	}

if( eta >2.5 ){EtaBin=2;
	if( pt <50  ) PtBin=1;
	if( 50<= pt && pt <80 ) PtBin=2;
	if( 80<= pt ) PtBin=3;
	}

float a,b;

long Selector=PtBin+10*RhoBin+100*EtaBin+ 1000*tag; // 10 is human readable

switch (Selector){
// tag eta rho pt
case 1111: a= .89; b= .13; break; 
case 1121: a= .94; b= -.03; break; 
case 1112: a= .91; b= -.08; break; 
case 1122: a= .96; b= -.28; break; 
case 1113: a= .865; b= -.05; break; 
case 1123: a= .94; b= .3; break; 

case 2111: a= .98; b= .04; break; 
case 2121: a= .99; b= -0.02; break; 
case 2112: a= .94; b= 0.14; break; 
case 2122: a= .995; b= 0.04; break; 
case 2113: a= .965; b= 0.08; break; 
case 2123: a= 1.030; b= -0.18; break; 

//FWD
case 1211: a= .89; b= -.42; break; 
case 1221: a= .795; b= -.38; break; 
case 1212: a= 1.075; b= -.66; break; 
case 1222: a= .986; b= -.236; break; 
case 1213: a= .95; b= -.30; break; 
case 1223: a= 1.03; b= -.32; break; 

case 2211: a= 1.12; b= 0.14; break; 
case 2221: a= 1.02; b= 0.020;  break;
case 2212: a= 1.10; b= 0.240;  break;
case 2222: a= .92; b= 0.2;  break;
case 2213: a= 1.00; b= 0.04;  break;
case 2223: a= 1.06; b= 0.1;  break;
default: return -1;
}

return function(value, a ,b);
}
