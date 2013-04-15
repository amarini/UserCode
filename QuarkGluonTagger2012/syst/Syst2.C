#include "TMath.h"
#include <string>

inline float function(float x0, float a ,float b,float min=0,float max=1)
{
//TMath::ATan( %f * TMath::Tan(TMath::Pi()*%s-TMath::Pi()/2.) + %f)/TMath::Pi() +0.5 
using namespace TMath;
//return ATan( a*Tan(Pi()*x-Pi()/2.0) + b)/Pi() +0.5;
//return (ATan( a*Tan( (Pi()*x)-(Pi()/2.0)) + b)/Pi()) + 0.5;

  //string var=Form("TMath::TanH( %f * TMath::ATanH(2*%s-1) + %f)/2+.5 ",par[0],varName.c_str(),par[1]);

float x=(x0-min)/(max-min); 
if(x<0)x=0;
if(x>1)x=1;

float x1= (TanH( a* ATanH(2*x-1)+b )/2+.5 ) ;

return x1*(max-min)+min;

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

float a=-999,b=-999,lmin=0,lmax=1;

int type=0;// 1 = PtBalance
long Selector = PtBin + 10*RhoBin + 100*EtaBin + 1000*tag ;//+ 10000*type; // 10 is human readable

switch (Selector){
//type Discr eta rho pt
case 3111: a=0.940; b=0.140; break;
case 3112: a=0.950; b=0.020; break;
case 3113: a=0.850; b=0.030; break;
case 3121: a=0.960; b=0.160; break;
case 3122: a=0.980; b=0.070; break;
case 3123: a=1.010; b=0.220; break;

case 3211: a=0.980; b=-0.04; break;
case 3212: a=0.920; b=-0.34; break;
case 3213: a=0.810; b=-0.16; break;
case 3221: a=0.910; b=-0.23; break;
case 3222: a=0.910; b=-0.10; break;
case 3223: a=0.900; b=-0.35; break;
//
case 2111: a=0.990; b=-0.06;lmin=0.075;lmax=0.876;break;
case 2112: a=0.980; b=0.040;lmin=0.065;lmax=0.900;break;
case 2113: a=0.960; b=0.060;lmin=0.060;lmax=0.917;break;
case 2121: a=1.010; b=-0.10;lmin=0.076;lmax=0.876;break;
case 2122: a=1.040; b=-0.08;lmin=0.066;lmax=0.901;break;
case 2123: a=1.060; b=-0.20;lmin=0.060;lmax=0.916;break;

case 2211: a=1.230; b=0.000;lmin=0.157;lmax=0.814;break;
case 2212: a=1.080; b=0.240;lmin=0.187;lmax=0.804;break;
case 2213: a=1.070; b=0.240;lmin=0.177;lmax=0.739;break;
case 2221: a=1.160; b=-0.08;lmin=0.170;lmax=0.817;break;
case 2222: a=1.080; b=0.200;lmin=0.188;lmax=0.864;break;
case 2223: a=1.040; b=0.120;lmin=0.177;lmax=0.817;break;

// tag eta rho pt
//MLP
//case 12111: a= .98; b= .03;lmin=0.075;lmax=0.876; break; 
//case 12112: a= .97; b= .12;lmin=0.065;lmax=0.90; break; 
//case 12113: a= .95; b= .08;lmin=0.06;lmax=0.917; break; 
//case 12121: a= .99; b= -.02;lmin=0.076;lmax=0.876; break; 
//case 12122: a= .94; b= .14;lmin=0.066;lmax=0.901; break; 
//case 12123: a= 1.030; b= -.090;lmin=0.06;lmax=0.916; break; 
//
////FWD - nPFCand -- for data
//
//case 12211: a=1.23; b=-0.000; lmin=.157; lmax=.814; break;
//case 12221: a=0.94; b=-0.090; lmin=.170; lmax=.817; break;
//
//case 12212: a=.910; b= 0.160; lmin=.187; lmax=.804; break;
//case 12222: a=.880; b= 0.180; lmin=.188; lmax=.864; break;
//
//case 12213: a=.840; b= 0.140; lmin=.177; lmax=.739; break;
//case 12223: a=1.130; b=0.020; lmin=.177; lmax=.817; break;
//
////HISTO
//// pt- 30-80 rho -0 -15
//case 13111: a=0.94; b=0.04; break;
//case 13121: a=0.95; b=0.06; break;
//case 13112: a=0.95; b=-0.14; break;
//case 13122: a=0.96; b=0.06 ; break;
//case 13113: a=0.87; b=-0.16; break;
//case 13123: a=.965; b=0.08; break;
////FWD - nPFCand -- for data
//case 13211: a=0.96; b=-.15; break;
//case 13221: a=.92; b=-.03; break;
//        //50-80
//case 13212: a=.88; b=-.35; break;
//case 13222: a=.95; b=.1; break;
//
//case 13213: a=.92; b=-.22; break;
//case 13223: a=.97; b=0.07; break;


default: return -1;
}

//if( (a<-990) && (b<-990)) printf("warning: a & b not set\n");
//printf("%.3f : %.3f : %.3f : %.3f\n",a,b,lmin,lmax);
return function(value, a ,b,lmin,lmax);
}
