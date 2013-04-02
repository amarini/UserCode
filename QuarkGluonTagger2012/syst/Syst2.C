#include "TMath.h"
#include <string>

inline float function(float x, float a ,float b)
{
//TMath::ATan( %f * TMath::Tan(TMath::Pi()*%s-TMath::Pi()/2.) + %f)/TMath::Pi() +0.5 
using namespace TMath;
//return ATan( a*Tan(Pi()*x-Pi()/2.0) + b)/Pi() +0.5;
//return (ATan( a*Tan( (Pi()*x)-(Pi()/2.0)) + b)/Pi()) + 0.5;

  //string var=Form("TMath::TanH( %f * TMath::ATanH(2*%s-1) + %f)/2+.5 ",par[0],varName.c_str(),par[1]);
return (TanH( a* ATanH(2*x-1)+b )/2+.5 ) ;

}

float Syst(const char *tagger,float pt, float rho, float eta,float value )//"QGL","MLP"
{
using namespace std;
int PtBin=-999;
int RhoBin=-999;
int EtaBin=-999;
int tag=-999;

if(string(tagger)==string("QGLFit")) 	tag=1;
if(string(tagger)==string("QGLMLP")) 	tag=2;
if(string(tagger)==string("MLP")) 	tag=2;
if(string(tagger)==string("QGLHisto")) 	tag=3;
if(string(tagger)==string("QGL")) 	tag=3;

if( rho <=15 ) RhoBin=1;
if( rho >15  ) RhoBin=2;

if( eta <=2.5) {EtaBin=1;
	if( pt <80  ) PtBin=1;
	if( 80<= pt && pt <120 ) PtBin=2;
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
//MLP
case 2111: a= .985; b= .02; break; 
case 2121: a= 1.00; b= .0; break; 
//HISTO
// pt- 30-80 rho -0 -15
case 3111: a=0.94; b=0.04; break;
case 3121: a=0.95; b=0.06; break;
case 3112: a=0.95; b=-0.14; break;
case 3122: a=0.96; b=0.06 ; break;
case 3113: a=0.87; b=-0.16; break;
case 3123: a=.965; b=0.08; break;
//FWD
case 3211: a=0.89; b=-.16; break;
case 3221: a=.83; b=-.22; break;
case 3212: a=.835; b=-.5; break;
case 3222: a=.765; b=-.18; break;
case 3213: a=.97; b=.04; break;
case 3223: a=.945; b=-.34; break;

default: return -1;
}

return function(value, a ,b);
}
