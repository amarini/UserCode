// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    Ntp1Finalizer_QG on QCD samples.
//
// ------------------------------------------------------------

#include <string>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"

#ifndef QGL
#define QGL

class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const char * nchargedFileName="ncharged_fitresults.root",const char * nneutralFileName="nneutral_fitresults.root",const char * PtDFileName="PtD_fitresults.root");
  virtual ~QGLikelihoodCalculator();

  float computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD );
//e' un po' una palla ma devo dichiararle qui
   static double gammadistr_(double* x, double* par)
	{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
	}
   static double functionPtD_(double * x ,double*par)
	{
        return TMath::Exp ( (x[0]-par[0])/par[1] *(
                                                -(x[0]-par[0])/par[1]+ TMath::Sqrt( TMath::Power( (x[0]-par[0])/par[1],2) +par[4] ) -par[2]
                                                )

        )//* par[3]
                * par[1]*(TMath::Sqrt(2*TMath::Pi())/2*par[3] +1./par[2] ); //normalizzazione a meta' (quasi 1 a parte per par3 che e' un'integrazione di una gauss)
	}
	   
 private:

  TFile *nchargedFile, *nneutralFile,*PtDFile;
  double x;
  double *nChargedPar,*nNeutralPar,*PtDPar;
  std::map<const char*,TGraph*> graph;

  float Interpolate(float pt, TGraph *g);
  //Fit functinos 
  

};

#endif
