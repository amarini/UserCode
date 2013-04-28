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

string Version="Syst3";

float Syst(const char *tagger,float pt, float rho, float eta,int pdgid,float value )//"QGL","MLP"
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
	if( pt <50  ) PtBin=1;
	if( 50 <= pt && pt <80  ) PtBin=2;
	if( 80<= pt && pt <120 ) PtBin=3;
	if( 120<= pt ) PtBin=4;
	}

if( eta >2.5 ){EtaBin=2;
	if( pt <50  ) PtBin=1;
	if( 50<= pt && pt <80 ) PtBin=2;
	if( 80<= pt ) PtBin=3;
	}

float a_q=-999,b_q=-999,a_g=-999,b_g=-999,lmin=0,lmax=1;

int type=0;// 1 = PtBalance
long Selector = PtBin + 10*RhoBin + 100*EtaBin + 1000*tag ;//+ 10000*type; // 10 is human readable

switch (Selector){
//type Discr eta rho pt
//QGLHisto: Pt=30_50 Rho=0_15 Eta=0_2
case 3111:a_q=0.910;b_q=0.040;a_g=1.000;b_g=0.120;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=30_50 Rho=15_40 Eta=0_2
case 3121:a_q=0.930;b_q=0.010;a_g=1.010;b_g=0.240;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=50_80 Rho=0_15 Eta=0_2
case 3112:a_q=0.900;b_q=0.060;a_g=0.970;b_g=-0.000;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=50_80 Rho=15_40 Eta=0_2
case 3122:a_q=0.920;b_q=0.020;a_g=1.010;b_g=0.140;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=80_120 Rho=0_15 Eta=0_2
case 3113:a_q=0.950;b_q=-0.180;a_g=0.870;b_g=-0.000;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=80_120 Rho=15_40 Eta=0_2
case 3123:a_q=0.890;b_q=0.040;a_g=0.910;b_g=-0.050;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=120_250 Rho=0_15 Eta=0_2
case 3114:a_q=0.850;b_q=-0.070;a_g=0.950;b_g=0.180;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=120_250 Rho=15_40 Eta=0_2
case 3124:a_q=0.990;b_q=-0.050;a_g=1.040;b_g=0.440;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=30_50 Rho=0_15 Eta=3_5
case 3211:a_q=1.010;b_q=-0.210;a_g=0.980;b_g=-0.010;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=30_50 Rho=15_40 Eta=3_5
case 3221:a_q=0.790;b_q=-0.140;a_g=1.090;b_g=0.090;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=50_80 Rho=0_15 Eta=3_5
case 3212:a_q=0.970;b_q=-0.270;a_g=0.980;b_g=-0.430;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=50_80 Rho=15_40 Eta=3_5
case 3222:a_q=0.940;b_q=-0.230;a_g=0.910;b_g=0.040;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=80_120 Rho=0_15 Eta=3_5
case 3213:a_q=0.980;b_q=-0.300;a_g=0.990;b_g=-0.010;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=80_120 Rho=15_40 Eta=3_5
case 3223:a_q=1.080;b_q=0.030;a_g=0.980;b_g=-0.450;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=120_250 Rho=0_15 Eta=3_5

default: return -1;
}

//if( (a<-990) && (b<-990)) printf("warning: a & b not set\n");
//printf("%.3f : %.3f : %.3f : %.3f\n",a,b,lmin,lmax);
if(pdgid==21)return function(value, a_g ,b_g,lmin,lmax);
else if ((fabs(pdgid)<4) && (pdgid!=0)) return function(value, a_q ,b_q,lmin,lmax);
else return value;
}
