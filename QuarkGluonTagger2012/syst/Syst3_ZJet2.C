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

string Version="Syst3.1";

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
long Selector = PtBin + 8*RhoBin + 8*8*EtaBin + 8*8*8*tag ;//+ 10000*type; // 10 is human readable

if( (Selector & 07000)==03000 ) switch (Selector){
//type Discr eta rho pt

//QGLHisto: Pt=30_50 Rho=0_15 Eta=0_2
case 03111:a_q=0.900;b_q=0.150;a_g=0.990;b_g=-0.000;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=30_50 Rho=15_40 Eta=0_2
case 03121:a_q=0.930;b_q=0.060;a_g=0.990;b_g=-0.000;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=50_80 Rho=0_15 Eta=0_2
case 03112:a_q=0.940;b_q=0.020;a_g=0.980;b_g=-0.000;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=50_80 Rho=15_40 Eta=0_2
case 03122:a_q=0.960;b_q=0.000;a_g=0.950;b_g=-0.030;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=80_120 Rho=0_15 Eta=0_2
case 03113:a_q=0.920;b_q=0.020;a_g=0.950;b_g=-0.050;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=80_120 Rho=15_40 Eta=0_2
case 03123:a_q=0.940;b_q=-0.000;a_g=0.930;b_g=-0.030;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=120_250 Rho=0_15 Eta=0_2
case 03114:a_q=0.940;b_q=-0.050;a_g=0.910;b_g=-0.080;lmin=0.000;lmax=1.000;break;
//QGLHisto: Pt=120_250 Rho=15_40 Eta=0_2
case 03124:a_q=0.940;b_q=-0.030;a_g=0.960;b_g=-0.010;lmin=0.000;lmax=1.000;break;


default: return -1;

} else if( (Selector&07000)==02000) switch(Selector){

//QGLMLP: Pt=30_50 Rho=0_15 Eta=0_2
case 02111:a_q=0.940;b_q=-0.020;a_g=0.920;b_g=-0.040;lmin=0.085;lmax=0.866;break;
//QGLMLP: Pt=30_50 Rho=15_40 Eta=0_2
case 02121:a_q=0.920;b_q=-0.030;a_g=0.910;b_g=-0.010;lmin=0.085;lmax=0.866;break;
//QGLMLP: Pt=50_80 Rho=0_15 Eta=0_2

//QGLMLP: Pt=80_120 Rho=0_15 Eta=0_2
case 02113:a_q=0.900;b_q=0.000;a_g=0.960;b_g=-0.080;lmin=0.065;lmax=0.910;break;
//QGLMLP: Pt=80_120 Rho=15_40 Eta=0_2
case 02123:a_q=0.890;b_q=-0.010;a_g=0.910;b_g=-0.060;lmin=0.065;lmax=0.902;break;
//QGLMLP: Pt=120_250 Rho=0_15 Eta=0_2

default: return -1;

}

//if( (a<-990) && (b<-990)) printf("warning: a & b not set\n");
//printf("%.3f : %.3f : %.3f : %.3f\n",a,b,lmin,lmax);
if(pdgid==21)return function(value, a_g ,b_g,lmin,lmax);
else if ((fabs(pdgid)<4) && (pdgid!=0)) return function(value, a_q ,b_q,lmin,lmax);
else return value;
}
