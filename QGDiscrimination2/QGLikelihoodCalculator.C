#include "QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TMath.h"
inline double QGLikelihoodCalculator::gammadistr_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}

inline double QGLikelihoodCalculator::functionPtD_(double * x ,double*par)
{
        return TMath::Exp ( (x[0]-par[0])/par[1] *(
                                                -(x[0]-par[0])/par[1]+ TMath::Sqrt( TMath::Power( (x[0]-par[0])/par[1],2) +par[4] ) -par[2]
                                                )

        ) //* par[3];
                * par[1]*(TMath::Sqrt(2*TMath::Pi())/2*par[3] +1./par[2] ); //normalizzazione a meta' (quasi 1 a parte per par3 che e' un'integrazione di una gauss)
}

//constructor
QGLikelihoodCalculator::QGLikelihoodCalculator(const char * nchargedFileName,const char * nneutralFileName,const char * PtDFileName)
{
nchargedFile=TFile::Open(nchargedFileName);
nneutralFile=TFile::Open(nneutralFileName);
PtDFile=TFile::Open(PtDFileName);
nChargedPar= new double[2];
nNeutralPar= new double[2];
PtDPar= new double[5];
}

float QGLikelihoodCalculator::computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD ) {
	double q,g;
	//setting parameters
	nChargedPar[0]=Interpolate(pt,(TGraph*)nchargedFile->Get("alpha_quark"));
	nChargedPar[1]=Interpolate(pt,(TGraph*)nchargedFile->Get("mean_quark"));
	nNeutralPar[0]=Interpolate(pt,(TGraph*)nneutralFile->Get("alpha_quark"));
	nNeutralPar[1]=Interpolate(pt,(TGraph*)nneutralFile->Get("mean_quark"));
	PtDPar[0]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD0_quark"));
	PtDPar[1]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD1_quark"));
	PtDPar[2]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD2_quark"));
	PtDPar[3]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD3_quark"));
	PtDPar[4]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD4_quark"));
	x=(double)nCharged;//it reads directly 
	q = gammadistr_(&x,nChargedPar);
	x=(double)nNeutral;
	q *= gammadistr_(&x,nNeutralPar);
	x=(double)ptD;
	q*= functionPtD_(&x,PtDPar);
	nChargedPar[0]=Interpolate(pt,(TGraph*)nchargedFile->Get("alpha_gluon"));
	nChargedPar[1]=Interpolate(pt,(TGraph*)nchargedFile->Get("mean_gluon"));
	nNeutralPar[0]=Interpolate(pt,(TGraph*)nneutralFile->Get("alpha_gluon"));
	nNeutralPar[1]=Interpolate(pt,(TGraph*)nneutralFile->Get("mean_gluon"));
	PtDPar[0]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD0_gluon"));
	PtDPar[1]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD1_gluon"));
	PtDPar[2]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD2_gluon"));
	PtDPar[3]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD3_gluon"));
	PtDPar[4]=Interpolate(pt,(TGraph*)PtDFile->Get("PtD4_gluon"));
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



