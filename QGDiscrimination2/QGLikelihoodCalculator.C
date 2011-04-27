#include "QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TMath.h"
//static double QGLikelihoodCalculator::gammadistr_(double* x, double* par)

//static double QGLikelihoodCalculator::functionPtD_(double * x ,double*par)

//constructor
using namespace std;
QGLikelihoodCalculator::QGLikelihoodCalculator(const char * nchargedFileName,const char * nneutralFileName,const char * PtDFileName)
{
nchargedFile=TFile::Open(nchargedFileName);
nneutralFile=TFile::Open(nneutralFileName);
PtDFile=TFile::Open(PtDFileName);
nChargedPar= new double[2];
nNeutralPar= new double[2];
PtDPar= new double[5];

//graph=new std::map<const char*,TGraph*>;
graph["alpha_quark_nc"]=(TGraph*)nchargedFile->Get("alpha_quark")->Clone("alpha_quark_nc");
graph["mean_quark_nc"]=(TGraph*)nchargedFile->Get("mean_quark")->Clone("mean_quark_nc");
graph["alpha_gluon_nc"]=(TGraph*)nchargedFile->Get("alpha_gluon")->Clone("alpha_gluon_nc");
graph["mean_gluon_nc"]=(TGraph*)nchargedFile->Get("mean_gluon")->Clone("mean_gluon_nc");

graph["alpha_quark_nn"]=(TGraph*)nneutralFile->Get("alpha_quark")->Clone("alpha_quark_nn");
graph["mean_quark_nn"]=(TGraph*)nneutralFile->Get("mean_quark")->Clone("mean_quark_nn");
graph["alpha_gluon_nn"]=(TGraph*)nneutralFile->Get("alpha_gluon")->Clone("alpha_gluon_nn");
graph["mean_gluon_nn"]=(TGraph*)nneutralFile->Get("mean_gluon")->Clone("mean_gluon_nn");

graph["PtD0_quark"]=(TGraph*)PtDFile->Get("PtD0_quark")->Clone("PtD0_quark");
graph["PtD1_quark"]=(TGraph*)PtDFile->Get("PtD1_quark")->Clone("PtD1_quark");
graph["PtD2_quark"]=(TGraph*)PtDFile->Get("PtD2_quark")->Clone("PtD2_quark");
graph["PtD3_quark"]=(TGraph*)PtDFile->Get("PtD3_quark")->Clone("PtD3_quark");
graph["PtD4_quark"]=(TGraph*)PtDFile->Get("PtD4_quark")->Clone("PtD4_quark");

graph["PtD0_gluon"]=(TGraph*)PtDFile->Get("PtD0_gluon")->Clone("PtD0_gluon");
graph["PtD1_gluon"]=(TGraph*)PtDFile->Get("PtD1_gluon")->Clone("PtD1_gluon");
graph["PtD2_gluon"]=(TGraph*)PtDFile->Get("PtD2_gluon")->Clone("PtD2_gluon");
graph["PtD3_gluon"]=(TGraph*)PtDFile->Get("PtD3_gluon")->Clone("PtD3_gluon");
graph["PtD4_gluon"]=(TGraph*)PtDFile->Get("PtD4_gluon")->Clone("PtD4_gluon");
}

QGLikelihoodCalculator::~QGLikelihoodCalculator()
{
nchargedFile->Close();
nneutralFile->Close();
PtDFile->Close();
delete nChargedPar;
delete nNeutralPar;
delete PtDPar;
//delete graph;
std::map<const char*,TGraph *>::iterator it;
for(it=graph.begin();it!=graph.end();it++)
{
	delete it->second;
}
}

float QGLikelihoodCalculator::computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD ) {
	double q,g;
	//setting parameters
	nChargedPar[0]=Interpolate(pt,graph["alpha_quark_nc"]);
	nChargedPar[1]=Interpolate(pt,graph["mean_quark_nc"]);
	nNeutralPar[0]=Interpolate(pt,graph["alpha_quark_nn"]);
	nNeutralPar[1]=Interpolate(pt,graph["mean_quark_nn"]);
	PtDPar[0]=Interpolate(pt,graph["PtD0_quark"]);
	PtDPar[1]=Interpolate(pt,graph["PtD1_quark"]);
	PtDPar[2]=Interpolate(pt,graph["PtD2_quark"]);
	PtDPar[3]=Interpolate(pt,graph["PtD3_quark"]);
	PtDPar[4]=Interpolate(pt,graph["PtD4_quark"]);
	x=(double)nCharged;//it reads directly 
	q = gammadistr_(&x,nChargedPar);
	x=(double)nNeutral;
	q *= gammadistr_(&x,nNeutralPar);
	x=(double)ptD;
	q*= functionPtD_(&x,PtDPar);
	nChargedPar[0]=Interpolate(pt,graph["alpha_gluon_nc"]);
	nChargedPar[1]=Interpolate(pt,graph["mean_gluon_nc"]);
	nNeutralPar[0]=Interpolate(pt,graph["alpha_gluon_nn"]);
	nNeutralPar[1]=Interpolate(pt,graph["mean_gluon_nn"]);
	PtDPar[0]=Interpolate(pt,graph["PtD0_gluon"]);
	PtDPar[1]=Interpolate(pt,graph["PtD1_gluon"]);
	PtDPar[2]=Interpolate(pt,graph["PtD2_gluon"]);
	PtDPar[3]=Interpolate(pt,graph["PtD3_gluon"]);
	PtDPar[4]=Interpolate(pt,graph["PtD4_gluon"]);
	x=(double)nCharged;//it reads directly 
	g = gammadistr_(&x,nChargedPar);
	x=(double)nNeutral;
	g *= gammadistr_(&x,nNeutralPar);
	x=(double)ptD;
	g*= functionPtD_(&x,PtDPar);
 return q/(q+g);
}


//Linear Interpolation
float QGLikelihoodCalculator::Interpolate(float pt, TGraph *g)
{
//suppongo che i bin siano ordinati - I'm supposing bin are sorted
Int_t Max=g->GetN()-1,Min=0;
double x0,y0,x1,y1;
//binary search
while(Max-Min>1)
	{
	g->GetPoint((Max-Min)/2+Min,x0,y0);
	if(pt<x0)Max=(Max-Min)/2+Min;
	else Min=(Max-Min)/2+Min;
	}
//ora ho in Min e Max i bin giusti
g->GetPoint(Min,x0,y0);
g->GetPoint(Max,x1,y1);
return (y0+ (y1-y0)/TMath::Log(x1/x0)*TMath::Log(pt/x0));

}



