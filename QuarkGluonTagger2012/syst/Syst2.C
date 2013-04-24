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

string Version="Syst2.1";

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
	if( pt <50  ) PtBin=1;
	if( 50< pt && pt <80  ) PtBin=2;
	if( 80<= pt && pt <120 ) PtBin=3;
	if( 120<= pt ) PtBin=4;
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

//QGLHisto: Pt=30_50 Rho=0_15 Eta=0_2
case 3111:a=0.930;b=0.040;lmin=0.000;lmax=1.000;break;//chi2=43.093; chi2_0=135.307
//QGLHisto: Pt=30_50 Rho=15_40 Eta=0_2
case 3121:a=0.960;b=0.080;lmin=0.000;lmax=1.000;break;//chi2=46.973; chi2_0=97.035
//QGLHisto: Pt=50_80 Rho=0_15 Eta=0_2
case 3112:a=0.930;b=0.050;lmin=0.000;lmax=1.000;break;//chi2=17.773; chi2_0=76.880
//QGLHisto: Pt=50_80 Rho=15_40 Eta=0_2
case 3122:a=0.950;b=0.060;lmin=0.000;lmax=1.000;break;//chi2=23.412; chi2_0=46.616
//QGLHisto: Pt=80_120 Rho=0_15 Eta=0_2
case 3113:a=0.920;b=-0.050;lmin=0.000;lmax=1.000;break;//chi2=17.317; chi2_0=68.765
//QGLHisto: Pt=80_120 Rho=15_40 Eta=0_2
case 3123:a=0.930;b=-0.010;lmin=0.000;lmax=1.000;break;//chi2=32.203; chi2_0=54.387
//QGLHisto: Pt=120_250 Rho=0_15 Eta=0_2
case 3114:a=0.820;b=-0.060;lmin=0.000;lmax=1.000;break;//chi2=25.969; chi2_0=134.824
//QGLHisto: Pt=120_250 Rho=15_40 Eta=0_2
case 3124:a=0.970;b=-0.040;lmin=0.000;lmax=1.000;break;//chi2=29.859; chi2_0=57.208
//QGLHisto: Pt=30_50 Rho=0_15 Eta=3_5
case 3211:a=0.980;b=-0.040;lmin=0.000;lmax=1.000;break;//chi2=33.541; chi2_0=46.062
//QGLHisto: Pt=30_50 Rho=15_40 Eta=3_5
case 3221:a=0.950;b=-0.100;lmin=0.000;lmax=1.000;break;//chi2=34.475; chi2_0=89.754
//QGLHisto: Pt=50_80 Rho=0_15 Eta=3_5
case 3212:a=0.910;b=-0.160;lmin=0.000;lmax=1.000;break;//chi2=48.519; chi2_0=109.996
//QGLHisto: Pt=50_80 Rho=15_40 Eta=3_5
case 3222:a=0.740;b=-0.200;lmin=0.000;lmax=1.000;break;//chi2=25.219; chi2_0=108.966
//QGLHisto: Pt=80_120 Rho=0_15 Eta=3_5
case 3213:a=0.810;b=-0.160;lmin=0.000;lmax=1.000;break;//chi2=64.986; chi2_0=150.560
//QGLHisto: Pt=80_120 Rho=15_40 Eta=3_5
case 3223:a=0.800;b=1.670;lmin=0.000;lmax=1.000;break;//chi2=38.111; chi2_0=83.738
//QGLHisto: Pt=120_250 Rho=0_15 Eta=3_5
case 3214:a=0.570;b=0.030;lmin=0.000;lmax=1.000;break;//chi2=38.924; chi2_0=70.384
//QGLHisto: Pt=120_250 Rho=15_40 Eta=3_5
case 3224:a=0.840;b=0.820;lmin=0.000;lmax=1.000;break;//chi2=5.228; chi2_0=10.585

//QGLMLP: Pt=30_50 Rho=0_15 Eta=0_2
case 2111:a=0.980;b=0.030;lmin=0.085;lmax=0.866;break;//chi2=32.503; chi2_0=52.354
//QGLMLP: Pt=30_50 Rho=15_40 Eta=0_2
case 2121:a=0.990;b=0.000;lmin=0.085;lmax=0.866;break;//chi2=29.040; chi2_0=36.774
//QGLMLP: Pt=50_80 Rho=0_15 Eta=0_2
case 2112:a=0.900;b=0.060;lmin=0.075;lmax=0.876;break;//chi2=209.055; chi2_0=445.726
//QGLMLP: Pt=50_80 Rho=15_40 Eta=0_2
case 2122:a=1.010;b=-0.040;lmin=0.076;lmax=0.876;break;//chi2=21.679; chi2_0=31.215
//QGLMLP: Pt=80_120 Rho=0_15 Eta=0_2
case 2113:a=0.950;b=0.070;lmin=0.065;lmax=0.900;break;//chi2=32.749; chi2_0=67.946
//QGLMLP: Pt=80_120 Rho=15_40 Eta=0_2
case 2123:a=1.000;b=0.010;lmin=0.066;lmax=0.901;break;//chi2=15.499; chi2_0=25.226
//QGLMLP: Pt=120_250 Rho=0_15 Eta=0_2
case 2114:a=0.920;b=0.150;lmin=0.060;lmax=0.917;break;//chi2=21.506; chi2_0=39.306
//QGLMLP: Pt=120_250 Rho=15_40 Eta=0_2
case 2124:a=1.010;b=-0.110;lmin=0.060;lmax=0.916;break;//chi2=18.821; chi2_0=37.154
//QGLMLP: Pt=30_50 Rho=0_15 Eta=3_5
case 2211:a=1.230;b=0.000;lmin=0.157;lmax=0.814;break;//chi2=29.988; chi2_0=149.568
//QGLMLP: Pt=30_50 Rho=15_40 Eta=3_5
case 2221:a=1.170;b=-0.080;lmin=0.170;lmax=0.831;break;//chi2=20.281; chi2_0=30.383
//QGLMLP: Pt=50_80 Rho=0_15 Eta=3_5
case 2212:a=1.040;b=0.340;lmin=0.187;lmax=0.804;break;//chi2=42.050; chi2_0=130.453
//QGLMLP: Pt=50_80 Rho=15_40 Eta=3_5
case 2222:a=1.080;b=0.200;lmin=0.188;lmax=0.864;break;//chi2=17.420; chi2_0=56.062
//QGLMLP: Pt=80_120 Rho=0_15 Eta=3_5
case 2213:a=1.140;b=0.320;lmin=0.177;lmax=0.739;break;//chi2=28.549; chi2_0=68.778
//QGLMLP: Pt=80_120 Rho=15_40 Eta=3_5
case 2223:a=0.000;b=1.500;lmin=0.177;lmax=0.817;break;//chi2=27.219; chi2_0=62.232
//QGLMLP: Pt=120_250 Rho=0_15 Eta=3_5

case 2214:a=0.820;b=0.530;lmin=0.180;lmax=0.707;break;//chi2=29.241; chi2_0=59.342
//QGLMLP: Pt=120_250 Rho=15_40 Eta=3_5
case 2224:a=1.010;b=1.510;lmin=0.189;lmax=0.667;break;//chi2=5.492; chi2_0=13.444
default: return -1;
}

//if( (a<-990) && (b<-990)) printf("warning: a & b not set\n");
//printf("%.3f : %.3f : %.3f : %.3f\n",a,b,lmin,lmax);
return function(value, a ,b,lmin,lmax);
}
